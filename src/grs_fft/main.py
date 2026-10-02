"""Ponto de entrada: `python -m grs_fft.main --radio vhf --iq-source tcp://grs-iq-rx:5556`."""

from __future__ import annotations

import argparse
import logging
import signal
import sys

from grs_fft.service import FftConfig, FftService


def parse_args(argv: list[str] | None = None) -> FftConfig:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radio", required=True,
                        help="Nome do rádio (vhf, uhf...): vira o tópico afc.<rádio> e fft.<rádio>")
    parser.add_argument("--iq-source", required=True, help="PUB de IQ do rádio (cf32_le, sem tópico)")
    parser.add_argument("--bind", default="tcp://*:5582", help="Onde publicar medidas e espectro")
    parser.add_argument("--sample-rate", type=float, default=240_000.0, help="Taxa do IQ, em S/s")
    parser.add_argument("--baud", type=float, default=1200.0,
                        help="Taxa do enlace: define a largura que um sinal precisa ter para "
                             "ser o satélite (1200 = beacon do FS-2, 4800 = dados)")
    parser.add_argument("--window-hz", type=float, default=8_000.0,
                        help="Até onde do centro procurar o sinal, em Hz")
    parser.add_argument("--min-snr-db", type=float, default=6.0,
                        help="SNR mínimo para uma medida valer")
    parser.add_argument("--spectrum-rate", type=float, default=5.0,
                        help="Quadros de espectro publicados por segundo, no máximo")
    args = parser.parse_args(argv)
    return FftConfig(radio=args.radio, iq_source=args.iq_source, bind=args.bind,
                     sample_rate_hz=args.sample_rate, baud=args.baud, window_hz=args.window_hz,
                     min_snr_db=args.min_snr_db, spectrum_rate_hz=args.spectrum_rate)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    try:
        config = parse_args(argv)
        service = FftService(config)
    except ValueError as error:
        print(f"[grs-fft] configuração inválida: {error}", file=sys.stderr)
        return 1

    signal.signal(signal.SIGTERM, lambda *_: service.stop())
    signal.signal(signal.SIGINT, lambda *_: service.stop())
    service.run()
    logging.getLogger(__name__).info("encerrado: %d quadros, %d rajadas medidas",
                                     service.frames, service.measurements)
    return 0


if __name__ == "__main__":
    sys.exit(main())
