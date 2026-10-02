"""Testes do bloco FFT.

O que se confere é o que a malha de ajuste vai usar: a medida cai no lugar do
sinal (a poucas dezenas de Hz), funciona longe do centro (onde o demodulador
já não decodifica), e NÃO aparece onde não há satélite — ruído, portadora
pura, sinal fora da janela. Uma medida falsa aqui vira o receptor sintonizado
no lugar errado.

O GFSK dos testes é gerado aqui, e não importado do simulador da estação:
medir com o mesmo código que gera prova menos do que parece.
"""

from __future__ import annotations

import json
import threading
import time

import numpy as np
import pytest
import zmq

from grs_fft.main import parse_args
from grs_fft.service import FftConfig, FftService
from grs_fft.spectrum import (
    BurstAggregator,
    Measurement,
    SpectrumAnalyzer,
    decimate_db,
    estimate_offset,
    gfsk_bandwidth_hz,
)

FS = 240_000


def gfsk(n_samples: int, baud: int, offset_hz: float = 0.0, h: float = 0.5,
         bt: float = 0.5, seed: int = 1) -> np.ndarray:
    """2GFSK com bits aleatórios, deslocado de `offset_hz`."""
    rng = np.random.default_rng(seed)
    sps = FS // baud
    bits = rng.integers(0, 2, n_samples // sps + 2) * 2 - 1
    nrz = np.repeat(bits, sps)[:n_samples].astype(float)
    sigma = np.sqrt(np.log(2)) / (2 * np.pi * bt) * sps
    taps = np.exp(-0.5 * (np.arange(-3 * sps, 3 * sps + 1) / sigma) ** 2)
    shaped = np.convolve(nrz, taps / taps.sum(), mode="same")
    deviation = h * baud / 2
    phase = 2 * np.pi * np.cumsum(deviation * shaped + offset_hz) / FS
    return np.exp(1j * phase).astype(np.complex64)


def noise(n_samples: int, power: float, seed: int = 2) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (np.sqrt(power / 2) * (rng.standard_normal(n_samples)
                                  + 1j * rng.standard_normal(n_samples))).astype(np.complex64)


def measure(samples: np.ndarray, baud: int, window_hz: float = 8000.0):
    analyzer = SpectrumAnalyzer(FS)
    frames = analyzer.feed(samples)
    return [estimate_offset(analyzer.freqs_hz, p, window_hz, gfsk_bandwidth_hz(baud))
            for p in frames]


# --- o espectro ----------------------------------------------------------------------


def test_tom_aparece_no_bin_certo_com_potencia_unitaria():
    analyzer = SpectrumAnalyzer(FS)
    t = np.arange(analyzer.frame_samples) / FS
    tone = np.exp(2j * np.pi * 15_000 * t).astype(np.complex64)

    power = analyzer.feed(tone)[0]

    assert analyzer.freqs_hz[np.argmax(power)] == pytest.approx(15_000, abs=analyzer.resolution_hz)
    assert np.max(power) == pytest.approx(1.0, rel=0.2)


def test_quadros_atravessam_lotes_sem_perder_amostras():
    analyzer = SpectrumAnalyzer(FS)
    data = gfsk(analyzer.frame_samples * 3, 1200)

    frames = []
    for chunk in np.array_split(data, 7):
        frames.extend(analyzer.feed(chunk))

    assert len(frames) == 3


# --- a medida --------------------------------------------------------------------------


@pytest.mark.parametrize("baud", [1200, 4800])
@pytest.mark.parametrize("offset", [0.0, 1200.0, -2500.0, 4700.0])
def test_mede_o_desvio_do_gfsk_longe_do_centro(baud, offset):
    """4,7 kHz: o pior oscilador declarado do TTC 2.0 em UHF. Bem além dos
    ~300 Hz em que o demodulador ainda decodifica."""
    n = SpectrumAnalyzer(FS).frame_samples * 4
    signal = gfsk(n, baud, offset) + noise(n, power=0.01)  # SNR ~20 dB na banda

    found = [m for m in measure(signal, baud) if m is not None]

    assert len(found) == 4
    assert np.mean([m.offset_hz for m in found]) == pytest.approx(offset, abs=40)


def test_mede_com_sinal_fraco():
    n = SpectrumAnalyzer(FS).frame_samples * 8
    # Potência do ruído na banda inteira: o SNR na banda do sinal fica ~10 dB.
    signal = gfsk(n, 1200, 800.0) + noise(n, power=13.0)

    found = [m for m in measure(signal, 1200) if m is not None]

    assert found, "nenhuma medida a ~10 dB"
    assert np.median([m.offset_hz for m in found]) == pytest.approx(800.0, abs=80)


def test_so_ruido_nao_vira_medida():
    n = SpectrumAnalyzer(FS).frame_samples * 20
    assert all(m is None for m in measure(noise(n, power=1.0), 1200))


def test_portadora_pura_nao_e_o_satelite():
    """Um espúrio ou interferência em CW, forte e dentro da janela, tem
    largura de um bin: travar nele mandaria o receptor para o lugar errado."""
    n = SpectrumAnalyzer(FS).frame_samples * 4
    t = np.arange(n) / FS
    carrier = np.exp(2j * np.pi * 2500 * t).astype(np.complex64)

    assert all(m is None for m in measure(carrier + noise(n, power=0.01), 1200))


def test_sinal_fora_da_janela_e_ignorado():
    n = SpectrumAnalyzer(FS).frame_samples * 4
    signal = gfsk(n, 1200, 20_000.0) + noise(n, power=0.01)

    assert all(m is None for m in measure(signal, 1200, window_hz=8000))


def test_sinal_largo_demais_e_recusado():
    """4800 baud medido por um bloco configurado para 1200: largura 4 vezes a
    esperada. Não é o enlace que este rádio deveria ouvir."""
    n = SpectrumAnalyzer(FS).frame_samples * 4
    signal = gfsk(n, 4800, 0.0) + noise(n, power=0.01)

    assert all(m is None for m in measure(signal, 1200))


def test_com_portadora_ao_lado_mede_o_gfsk():
    n = SpectrumAnalyzer(FS).frame_samples * 4
    t = np.arange(n) / FS
    carrier = 0.3 * np.exp(2j * np.pi * (-5000) * t).astype(np.complex64)
    signal = gfsk(n, 1200, 1500.0) + carrier + noise(n, power=0.01)

    found = [m for m in measure(signal, 1200) if m is not None]

    assert found and np.mean([m.offset_hz for m in found]) == pytest.approx(1500, abs=60)


# --- rajadas -----------------------------------------------------------------------------


def m(offset, power=1.0):
    return Measurement(offset_hz=offset, snr_db=20.0, bandwidth_hz=1800.0, excess_power=power)


def test_uma_medida_por_rajada_ponderada_pela_potencia():
    burst = BurstAggregator()

    assert burst.push(m(1000, 1.0)) is None
    assert burst.push(m(1200, 3.0)) is None
    closed = burst.push(None)

    assert closed.frames == 2
    assert closed.offset_hz == pytest.approx(1150.0)
    assert burst.push(None) is None


def test_rajada_longa_e_fechada_em_pedacos():
    burst = BurstAggregator(max_frames=3)

    results = [burst.push(m(100)) for _ in range(3)]

    assert results[:2] == [None, None] and results[2].frames == 3


def test_decimate_db_reduz_o_espectro():
    power = np.ones(4096)
    out = decimate_db(power, 512)

    assert len(out) == 512 and out[0] == pytest.approx(0.0)


# --- serviço por ZMQ ----------------------------------------------------------------------


def test_servico_publica_medida_por_rajada_e_espectro():
    context = zmq.Context()
    iq = context.socket(zmq.PUB)
    iq_port = iq.bind_to_random_port("tcp://127.0.0.1")
    service = FftService(FftConfig(radio="vhf", iq_source=f"tcp://127.0.0.1:{iq_port}",
                                   bind="tcp://127.0.0.1:0", spectrum_rate_hz=100.0), context)
    out = context.socket(zmq.SUB)
    out.setsockopt(zmq.SUBSCRIBE, b"afc.vhf")
    out.setsockopt(zmq.SUBSCRIBE, b"fft.vhf")
    out.setsockopt(zmq.RCVTIMEO, 5000)
    out.connect(service.endpoint)
    thread = threading.Thread(target=service.run, daemon=True)
    thread.start()
    time.sleep(0.5)

    frame = SpectrumAnalyzer(FS).frame_samples
    burst = gfsk(frame * 3, 1200, 900.0) + noise(frame * 3, power=0.01)
    silence = noise(frame * 2, power=0.01, seed=5)
    try:
        for block in np.array_split(np.concatenate((burst, silence)), 10):
            iq.send(block.tobytes())
            time.sleep(0.02)

        seen = {}
        deadline = time.time() + 5
        while time.time() < deadline and "afc.vhf" not in seen:
            parts = out.recv_multipart()
            seen[parts[0].decode()] = parts
    finally:
        service.stop()
        thread.join(timeout=3)
        out.close()
        iq.close()
        context.term()

    afc = json.loads(seen["afc.vhf"][1])
    assert afc["radio"] == "vhf" and afc["offset_hz"] == pytest.approx(900, abs=40)
    header = json.loads(seen["fft.vhf"][1])
    assert header["bins"] == 512
    assert len(np.frombuffer(seen["fft.vhf"][2], dtype=np.float32)) == 512


@pytest.mark.parametrize("changes, message", [
    ({"radio": "v h f"}, "rádio"),
    ({"window_hz": 200_000}, "janela"),
    ({"baud": 0}, "baud"),
])
def test_configuracao_invalida_e_recusada(changes, message):
    config = FftConfig(**{"radio": "vhf", "iq_source": "tcp://x:1", **changes})

    with pytest.raises(ValueError, match=message):
        config.validate()


def test_linha_de_comando():
    config = parse_args(["--radio", "uhf", "--iq-source", "tcp://grs-iq-rx-uhf:5556",
                         "--baud", "4800"])

    assert (config.radio, config.baud, config.window_hz) == ("uhf", 4800, 8000)



def test_dados_desequilibrados_enviesam_a_medida():
    """A premissa do centro de massa: tantos 0 quanto 1. Com 37,5% de uns (o
    payload `00 01 02 ... 3F` sem embaralhar), o tom do 0 pesa mais e a medida
    escorrega para o lado dele. O NGHam embaralha (CCSDS) e não tem isso; um
    enlace sem scrambler teria."""
    rng = np.random.default_rng(3)
    n = SpectrumAnalyzer(FS).frame_samples * 4
    sps = FS // 1200
    bits = np.where(rng.random(n // sps + 2) < 0.375, 1, -1)
    nrz = np.repeat(bits, sps)[:n].astype(float)
    sigma = np.sqrt(np.log(2)) / (2 * np.pi * 0.5) * sps
    taps = np.exp(-0.5 * (np.arange(-3 * sps, 3 * sps + 1) / sigma) ** 2)
    shaped = np.convolve(nrz, taps / taps.sum(), mode="same")
    signal = np.exp(1j * 2 * np.pi * np.cumsum(300 * shaped) / FS).astype(np.complex64)

    found = [m for m in measure(signal + noise(n, power=0.01), 1200) if m is not None]

    bias = np.mean([m.offset_hz for m in found])
    assert bias == pytest.approx((2 * 0.375 - 1) * 300, abs=30)
