# Contexto do projeto

Bloco FFT da estação terrestre SpaceLab ("FFT x3" no diagrama), um por rádio.
Mede onde o sinal do satélite está em relação à sintonia; o Station Manager
fecha a malha de ajuste fino (AFC) com essa medida. Ver o README.

## Desenho

```
spectrum.py   DSP puro: espectro de Welch, medida do desvio, rajadas
service.py    ZMQ: SUB do IQ, PUB de afc.<rádio> e fft.<rádio>
main.py       linha de comando
```

## Decisões, e por quê

**O IQ não sai do Station Server.** São ~1,9 MB/s por rádio. Este bloco mora
ao lado do receptor exatamente para transformar isso em bytes por rajada
(`afc`) e dezenas de KB/s (`fft`) — o que atravessa a rede até o Station
Manager.

**Quem integra é o Station Manager, não este bloco.** Ele é a única
autoridade sobre a sintonia: sabe se há passagem, de que satélite, em que
rádio, e zera a malha a cada passagem. Este bloco só mede e diz o que viu.

**Centro de massa do espectro, e não o demodulador.** O rastreador de DC do
demodulador também estima o desvio, mas só dentro de ~300 Hz — com o
oscilador do satélite a 1,5 kHz, o demodulador não decodifica e não mede. A
FFT acha o sinal em qualquer ponto da janela.

**Uma medida por rajada, não por quadro.** Um quadro (68 ms) pega a rajada
pela metade; cada pedaço integrado como medida inteira faria a malha andar
demais.

## Armadilhas

- **PREMISSA: dados equilibrados.** O centro de massa cai na portadora porque
  os dois tons do GFSK pesam igual. O NGHam do FS-2 embaralha o codeword com
  a sequência CCSDS, e isso garante. Sem scrambler, uma fração p de uns
  enviesa a medida em (2p − 1)·desvio: o payload `00 01 02 ... 3F` do
  simulador, sem embaralhar, deu −50 Hz a 1200 baud — e por isso o simulador
  passou a embaralhar. O teste `test_dados_desequilibrados_enviesam_a_medida`
  mostra o tamanho.
- **A 4800 baud a medida por rajada varia mais** (±30–60 Hz no simulador): a
  rajada é curta (0,17 s, 2–3 quadros) e cada uma é cortada de um jeito. A
  média entre rajadas, no Station Manager, cuida disso.
- **A `--baud` define o que é "o satélite".** Configurada errada, a largura
  esperada erra e toda medida é recusada (ou pior, uma interferência passa).

## Convenções

- Comentários e mensagens de commit em português; código em inglês.
- Testes conferem propriedades do sinal (onde a medida cai, o que é
  recusado), não que o processo roda. O GFSK dos testes é gerado aqui, e não
  importado do simulador da estação.
