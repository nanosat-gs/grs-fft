"""O bloco FFT de um rádio: assina o IQ, publica medidas e espectro.

    SUB  <fonte de IQ>    cf32_le, um lote por mensagem, sem tópico (o
                          contrato da :5556 do grs-iq-rx)
    PUB  <bind>           [b"afc.<rádio>", JSON]
                              uma por RAJADA com sinal: o desvio do sinal em
                              relação ao centro da sintonia (o resíduo, depois
                              do Doppler previsto), SNR, largura, quadros
                          [b"fft.<rádio>", JSON, float32[]]
                              o espectro reduzido, em dB, a no máximo
                              `spectrum_rate_hz` quadros por segundo

O IQ não sai desta máquina: são ~1,9 MB/s por rádio. O que atravessa a rede
até o Station Manager são estas duas mensagens — bytes por rajada e algumas
dezenas de KB/s de espectro.

Quem decide o que fazer com a medida é o Station Manager (a malha de ajuste
fica lá): ele sabe se há passagem, de que satélite, em que rádio. Este bloco
só mede e diz o que viu.
"""

from __future__ import annotations

import json
import logging
import threading
import time
from dataclasses import dataclass

import numpy as np
import zmq

from grs_fft.spectrum import (
    BurstAggregator,
    SpectrumAnalyzer,
    decimate_db,
    estimate_offset,
    gfsk_bandwidth_hz,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class FftConfig:
    radio: str
    iq_source: str
    bind: str = "tcp://*:5582"
    sample_rate_hz: float = 240_000.0
    baud: float = 1200.0
    modulation_index: float = 0.5
    # Quanto longe do centro procurar o sinal. Precisa cobrir o erro do
    # oscilador do satélite (±10 ppm declarados no TTC 2.0: ±1,5 kHz em VHF,
    # ±4,7 kHz em UHF) mais o que sobra do Doppler previsto.
    window_hz: float = 8_000.0
    min_snr_db: float = 6.0
    fft_size: int = 4096
    averages: int = 4
    spectrum_bins: int = 512
    spectrum_rate_hz: float = 5.0

    def validate(self) -> None:
        if not self.radio or not self.radio.replace("-", "").replace("_", "").isalnum():
            raise ValueError(f"nome de rádio inválido: {self.radio!r}")
        if self.window_hz <= 0 or self.window_hz >= self.sample_rate_hz / 2:
            raise ValueError("window_hz precisa caber na banda (0 < janela < taxa/2)")
        if self.baud <= 0:
            raise ValueError("baud precisa ser positivo")
        if self.fft_size % self.spectrum_bins:
            raise ValueError("spectrum_bins precisa dividir fft_size")


class FftService:
    def __init__(self, config: FftConfig, context: zmq.Context | None = None) -> None:
        config.validate()
        self.config = config
        self._context = context or zmq.Context.instance()
        self._analyzer = SpectrumAnalyzer(config.sample_rate_hz, config.fft_size, config.averages)
        self._bursts = BurstAggregator()
        self._expected_bw = gfsk_bandwidth_hz(config.baud, config.modulation_index)
        self._afc_topic = f"afc.{config.radio}".encode()
        self._fft_topic = f"fft.{config.radio}".encode()
        self._last_spectrum = 0.0
        self._stop = threading.Event()

        self.frames = 0
        self.measurements = 0

        self._sub = self._context.socket(zmq.SUB)
        self._sub.setsockopt(zmq.SUBSCRIBE, b"")
        self._sub.setsockopt(zmq.RCVTIMEO, 300)
        self._sub.setsockopt(zmq.LINGER, 0)
        self._sub.connect(config.iq_source)

        self._pub = self._context.socket(zmq.PUB)
        self._pub.setsockopt(zmq.LINGER, 0)
        self._pub.bind(config.bind)

    @property
    def endpoint(self) -> str:
        return self._pub.getsockopt_string(zmq.LAST_ENDPOINT)

    def stop(self) -> None:
        self._stop.set()

    def process(self, samples: np.ndarray) -> None:
        """Um lote de IQ: fecha quadros, mede, publica. Público para teste."""
        for power in self._analyzer.feed(samples):
            self.frames += 1
            measurement = estimate_offset(
                self._analyzer.freqs_hz, power, self.config.window_hz,
                self._expected_bw, self.config.min_snr_db,
            )
            self._publish_burst(self._bursts.push(measurement))
            self._publish_spectrum(power)

    def flush(self) -> None:
        """Sem IQ chegando: fecha a rajada em curso, se houver."""
        self._publish_burst(self._bursts.push(None))

    def run(self) -> None:
        logger.info("FFT do rádio %s: %s -> %s | janela ±%.0f Hz, GFSK de %.0f Hz "
                    "(%.0f baud), quadro de %.0f ms, resolução %.1f Hz",
                    self.config.radio, self.config.iq_source, self.config.bind,
                    self.config.window_hz, self._expected_bw, self.config.baud,
                    self._analyzer.frame_seconds * 1000, self._analyzer.resolution_hz)
        try:
            while not self._stop.is_set():
                try:
                    message = self._sub.recv()
                except zmq.Again:
                    self.flush()
                    continue
                except zmq.ZMQError:
                    if self._stop.is_set():
                        break
                    raise
                if len(message) % 8:
                    logger.warning("lote de IQ com %d bytes (não é cf32): descartado", len(message))
                    continue
                self.process(np.frombuffer(message, dtype=np.complex64))
        finally:
            self._sub.close()
            self._pub.close()

    def _publish_burst(self, burst) -> None:
        if burst is None:
            return
        self.measurements += 1
        payload = {"radio": self.config.radio, "t": time.time(), **burst.to_dict()}
        self._pub.send_multipart([self._afc_topic, json.dumps(payload).encode()], zmq.NOBLOCK)
        logger.info("rajada no %s: %+.0f Hz do centro, SNR %.1f dB, largura %.0f Hz, %d quadros",
                    self.config.radio, burst.offset_hz, burst.snr_db, burst.bandwidth_hz,
                    burst.frames)

    def _publish_spectrum(self, power: np.ndarray) -> None:
        now = time.monotonic()
        if now - self._last_spectrum < 1.0 / self.config.spectrum_rate_hz:
            return
        self._last_spectrum = now
        header = {"radio": self.config.radio, "t": time.time(),
                  "sample_rate_hz": self.config.sample_rate_hz,
                  "bins": self.config.spectrum_bins}
        values = np.asarray(decimate_db(power, self.config.spectrum_bins), dtype=np.float32)
        self._pub.send_multipart([self._fft_topic, json.dumps(header).encode(), values.tobytes()],
                                 zmq.NOBLOCK)
