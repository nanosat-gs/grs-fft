"""Bloco FFT da estação terrestre SpaceLab ("FFT x3" no diagrama).

Um por rádio, ao lado do receptor: assina o IQ, calcula o espectro em volta
da sintonia e mede a que distância do centro o sinal do satélite está. A
medida vai ao Station Manager, que fecha a malha de ajuste fino da sintonia;
o espectro vai ao Station Manager, que o repassa a quem desenha.
"""
