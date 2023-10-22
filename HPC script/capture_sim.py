from mps_capture_simulation import Sistema
import imageio
from tqdm import tqdm
import numpy as np

ENERGIA_INTERNA = 1000000
VOLUME = 360000
NUMERO_DE_PARTICULAS_1 = 100
NUMERO_DE_PARTICULAS_2 = 30
RAIO_PARTICULA_1 = 2
RAIO_PARTICULA_2 = 5

Sistema_particulas = Sistema(ENERGIA_INTERNA, VOLUME, NUMERO_DE_PARTICULAS_1, NUMERO_DE_PARTICULAS_2, RAIO_PARTICULA_1, RAIO_PARTICULA_2)

N = 100

for i in tqdm(range(N)):
    Sistema_particulas.main(0.01, 10, i, N)

frames = []
for i in range(N):
    image = imageio.imread(f'../GIF - mps capture simulation/Imagem {i}.jpg')
    frames.append(image)

imageio.mimsave('./mps capture simulation.gif', # output gif
               frames,          # array of input frames
               duration = 0.1)         # optional: frames per second

concs = np.array([Sistema_particulas.concentracoes['A'], Sistema_particulas.concentracoes['B'], Sistema_particulas.concentracoes['C']])
np.save('./capturas.npy', concs)
