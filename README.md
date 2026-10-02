# GRS FFT

Bloco FFT da estação terrestre SpaceLab (o "FFT x3" do diagrama da estação):
**um por rádio**, ao lado do receptor. Assina o IQ do rádio, calcula o espectro
em volta da sintonia e mede, a cada rajada do satélite, a que distância do
centro ela chegou.

```
                     ┌───────────── Station Server ─────────────┐     ┌── Control Server ──┐
  SDR ─▶ IQ Receiver ─┬─▶ demodulador ─▶ ...                      │     │                    │
                      └─▶ [ grs-fft ] ── afc.<rádio> (medida) ────┼────▶│  Station Manager   │
                            (IQ não sai daqui)  fft.<rádio> ──────┼────▶│  (malha + repasse) │
                      ◀── tune ── sintetizador ◀── offset.<rádio> ┼─────│                    │
                     └───────────────────────────────────────────┘     └────────────────────┘
```

## Por que existe

O Doppler que a estação corrige é **previsto** pelo TLE. O que a previsão não
vê é o erro do oscilador do satélite (o rádio do FS-2 declara cristal de ±10
ppm: até ±1,5 kHz em 145,9 MHz e ±4,7 kHz em 468,4 MHz) e o erro do próprio
TLE. O demodulador só tolera ~300 Hz de erro de sintonia: com o oscilador fora
assim, ele não decodifica nada — e não tem como medir nada. O espectro, sim:
olha uma janela larga e acha o sinal onde quer que ele esteja.

Medido na estação, com o simulador imitando a ISS e erro de oscilador de +1200
Hz (beacon VHF) e −3000 Hz (dados UHF): só com o Doppler previsto, 0 pacotes
nas duas cadeias; com o bloco FFT fechando a malha, o ajuste aprendeu +1194 a
+1203 Hz e −2971 Hz, os sinais voltaram para a poucas dezenas de hertz do
centro e os pacotes voltaram.

## O que publica

    SUB  <fonte de IQ>    cf32_le, um lote por mensagem, sem tópico
    PUB  <bind> (:5582)   [b"afc.<rádio>", JSON]
                              uma por rajada: offset_hz (o resíduo, depois do
                              Doppler previsto), snr_db, bandwidth_hz, frames
                          [b"fft.<rádio>", JSON, float32[512]]
                              o espectro em dB, até 5 quadros por segundo

Quem integra a medida é o **Station Manager** (`--fft-sources`), não este
bloco: ele sabe se há passagem, de que satélite, em que rádio, e é a única
autoridade sobre a sintonia. Ele também repassa o `fft.*` num endereço só
(`--spectrum-bind`), para o Spectrum Monitor.

## Usando

```bash
pip install -e ".[dev]"
pytest

python -m grs_fft.main --radio vhf --iq-source tcp://grs-iq-rx:5556 --baud 1200
python -m grs_fft.main --radio uhf --iq-source tcp://grs-iq-rx-uhf:5556 --baud 4800
```

Na estação (`nanosat-gs/grs-station`) sobem os serviços `grs-fft` e
`grs-fft-uhf`, nos profiles `rx` e `rxsim`.

## O que ele recusa, e por quê

Travar numa coisa que não é o satélite seria pior que não corrigir. Uma medida
só vale se o sinal:

- está dentro da janela de busca (`--window-hz`, 8 kHz por padrão);
- passa do SNR mínimo (`--min-snr-db`, 6 dB);
- tem a largura de um 2GFSK naquela taxa (Carson: `baud · (1 + h)`, aceita de
  0,3 a 2,5 vezes) — uma portadora pura (interferência, espúrio) é estreita
  demais; outro enlace, ou ruído de banda, largo demais.
