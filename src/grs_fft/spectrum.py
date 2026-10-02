"""O espectro em volta da sintonia, e onde o sinal do satélite está nele.

Duas contas, as duas sobre o mesmo periodograma:

- `SpectrumAnalyzer` acumula IQ e devolve o espectro de potência médio de um
  quadro (Welch: segmentos com janela de Hann, média das potências). Um FFT
  único piscaria: um espinho de ruído venceria um sinal fraco de verdade.
- `estimate_offset` acha a rajada do satélite num quadro e mede a que
  distância do centro ela está.

O desvio medido é o RESÍDUO: o receptor já está sintonizado em nominal +
Doppler previsto (+ o que a malha já corrigiu), então o que sobra é o erro do
oscilador do satélite e do TLE. É isso que o Station Manager integra.

Por que o centro de massa do espectro: num 2GFSK com dados equilibrados os
dois tons aparecem com a mesma potência, simétricos em volta da portadora, e o
centro de massa cai nela. PREMISSA: dados equilibrados. O NGHam do FS-2
embaralha o codeword com a sequência CCSDS, e isso os equilibra; num enlace
sem scrambler, uma fração p de uns enviesa a medida em (2p - 1)·desvio
(37,5% de uns a 1200 baud: −75 Hz). O teste
`test_dados_desequilibrados_enviesam_a_medida` mostra o tamanho. Funciona com o sinal a quilohertz do centro — onde
o demodulador já não decodifica nada e, portanto, não tem o que medir.

O que NÃO é sinal do satélite, e é recusado:
- ruído: abaixo do SNR mínimo;
- portadora pura (interferência, espúrio do próprio rádio): largura muito
  menor que a de um GFSK naquela taxa;
- coisa larga demais (outro serviço, ruído de banda): largura muito maior;
- qualquer coisa fora da janela de busca em volta do centro.
Travar numa portadora que não é o satélite seria pior que não corrigir.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class SpectrumAnalyzer:
    """Acumula IQ e entrega um espectro médio por quadro."""

    def __init__(self, sample_rate_hz: float, fft_size: int = 4096, averages: int = 4) -> None:
        if fft_size < 64 or fft_size & (fft_size - 1):
            raise ValueError("fft_size precisa ser potência de 2, >= 64")
        if averages < 1:
            raise ValueError("averages precisa ser >= 1")
        self.sample_rate_hz = float(sample_rate_hz)
        self.fft_size = fft_size
        self.averages = averages
        self._window = np.hanning(fft_size).astype(np.float32)
        # Normalização: potência de um tom de amplitude 1 não depende do tamanho.
        self._scale = 1.0 / float(np.sum(self._window) ** 2)
        self._pending = np.zeros(0, dtype=np.complex64)
        self.freqs_hz = np.fft.fftshift(np.fft.fftfreq(fft_size, 1.0 / self.sample_rate_hz))

    @property
    def frame_samples(self) -> int:
        return self.fft_size * self.averages

    @property
    def frame_seconds(self) -> float:
        return self.frame_samples / self.sample_rate_hz

    @property
    def resolution_hz(self) -> float:
        return self.sample_rate_hz / self.fft_size

    def feed(self, samples: np.ndarray) -> list[np.ndarray]:
        """Devolve o espectro (potência linear, de -fs/2 a +fs/2) de cada
        quadro completo que estas amostras fecharam."""
        data = np.concatenate((self._pending, samples.astype(np.complex64, copy=False)))
        frames = []
        count = len(data) // self.frame_samples
        for k in range(count):
            chunk = data[k * self.frame_samples:(k + 1) * self.frame_samples]
            segments = chunk.reshape(self.averages, self.fft_size) * self._window
            power = np.mean(np.abs(np.fft.fft(segments, axis=1)) ** 2, axis=0) * self._scale
            frames.append(np.fft.fftshift(power))
        self._pending = data[count * self.frame_samples:]
        return frames


def gfsk_bandwidth_hz(baud: float, modulation_index: float = 0.5) -> float:
    """Largura de Carson de um 2GFSK: 2·desvio + taxa = taxa·(1 + h).

    1200 baud (beacon do FS-2): 1800 Hz. 4800 baud (dados): 7200 Hz.
    """
    return baud * (1.0 + modulation_index)


@dataclass(frozen=True)
class Measurement:
    """Onde o sinal está num quadro, e quão confiável é."""

    offset_hz: float
    snr_db: float
    bandwidth_hz: float
    # Potência acima do ruído: o peso desta medida na média da rajada.
    excess_power: float


def estimate_offset(
    freqs_hz: np.ndarray,
    power: np.ndarray,
    window_hz: float,
    expected_bandwidth_hz: float,
    min_snr_db: float = 6.0,
    bandwidth_range: tuple[float, float] = (0.3, 2.5),
) -> Measurement | None:
    """O sinal do satélite neste quadro, ou None se não há um.

    A região do sinal é o trecho contíguo em volta do pico (dentro da janela)
    acima do limiar de SNR, depois de suavizar o espectro na escala de meia
    largura esperada — sem a suavização, os dois tons do GFSK com o vale entre
    eles virariam duas regiões, e a de cada tom seria estreita como uma
    portadora.
    """
    resolution = float(freqs_hz[1] - freqs_hz[0])
    # Mediana do espectro inteiro: o sinal ocupa uma fração pequena da banda
    # (7,2 kHz em 240 kHz no pior caso), então ela é o ruído.
    noise = float(np.median(power))
    if noise <= 0:
        return None

    inside = np.abs(freqs_hz) <= window_hz
    smooth_bins = max(1, int(round(expected_bandwidth_hz / 2 / resolution)))
    smoothed = np.convolve(power, np.ones(smooth_bins) / smooth_bins, mode="same")

    threshold = noise * 10.0 ** (min_snr_db / 10.0)
    candidates = np.flatnonzero(inside & (smoothed > threshold))
    if candidates.size == 0:
        return None

    peak = candidates[np.argmax(smoothed[candidates])]
    low = high = peak
    while low - 1 >= 0 and inside[low - 1] and smoothed[low - 1] > threshold:
        low -= 1
    while high + 1 < len(power) and inside[high + 1] and smoothed[high + 1] > threshold:
        high += 1

    # A suavização alarga a região em ~meia janela de cada lado; desconta.
    bandwidth = max(0.0, (high - low + 1) * resolution - (smooth_bins - 1) * resolution)
    low_bw, high_bw = bandwidth_range
    if not low_bw * expected_bandwidth_hz <= bandwidth <= high_bw * expected_bandwidth_hz:
        return None

    region = slice(low, high + 1)
    excess = np.clip(power[region] - noise, 0.0, None)
    total = float(np.sum(excess))
    if total <= 0:
        return None

    offset = float(np.sum(freqs_hz[region] * excess) / total)
    snr = 10.0 * np.log10(float(np.mean(power[region])) / noise)
    return Measurement(offset_hz=offset, snr_db=snr, bandwidth_hz=bandwidth, excess_power=total)


@dataclass(frozen=True)
class BurstMeasurement:
    """Uma rajada inteira: a média ponderada das medidas dos seus quadros."""

    offset_hz: float
    snr_db: float
    bandwidth_hz: float
    frames: int

    def to_dict(self) -> dict:
        return {"offset_hz": round(self.offset_hz, 1), "snr_db": round(self.snr_db, 1),
                "bandwidth_hz": round(self.bandwidth_hz, 1), "frames": self.frames}


class BurstAggregator:
    """Junta os quadros consecutivos com sinal numa medida por rajada.

    Uma por rajada, e não uma por quadro: um quadro pega a rajada pela metade
    (começo ou fim), e o Station Manager integraria cada pedaço como se fosse
    uma medida inteira. `max_frames` fecha uma rajada longa demais para a
    malha não ficar esperando um sinal contínuo terminar.
    """

    def __init__(self, max_frames: int = 30) -> None:
        self._max_frames = max_frames
        self._items: list[Measurement] = []

    def push(self, measurement: Measurement | None) -> BurstMeasurement | None:
        if measurement is not None:
            self._items.append(measurement)
            if len(self._items) < self._max_frames:
                return None
        return self._close()

    def _close(self) -> BurstMeasurement | None:
        items, self._items = self._items, []
        if not items:
            return None
        weights = np.array([m.excess_power for m in items])
        offsets = np.array([m.offset_hz for m in items])
        return BurstMeasurement(
            offset_hz=float(np.sum(offsets * weights) / np.sum(weights)),
            snr_db=float(max(m.snr_db for m in items)),
            bandwidth_hz=float(np.median([m.bandwidth_hz for m in items])),
            frames=len(items),
        )


def decimate_db(power: np.ndarray, bins: int) -> list[float]:
    """O espectro em `bins` faixas, em dB — o tamanho que vale a pena mandar
    pela rede para alguém desenhar."""
    if len(power) % bins:
        raise ValueError("bins precisa dividir o tamanho do FFT")
    grouped = power.reshape(bins, -1).mean(axis=1)
    return np.round(10.0 * np.log10(grouped + 1e-20), 1).tolist()
