# simulador-rlc-rungekutta
Simulador interativo de um circuito **RLC série**, resolvido numericamente via **método de Runge–Kutta de 4ª ordem**, com visualização dinâmica da corrente \( I(t) \) ao longo do tempo. 

## Funcionalidades

- Ajuste interativo dos parâmetros do circuito:
  - Resistência \( R \) (Ω)
  - Indutância \( L \) (H)
  - Capacitância \( C \) (F)
- Seleção do tipo de fonte:
  - Corrente contínua (DC)
  - Corrente alternada (AC senoidal)
- Controle de:
  - Tensão da fonte DC
  - Amplitude e frequência da fonte AC
- Cálculo e exibição de:
  - Frequência de ressonância
  - Impedância total do circuito
  - Ângulo de fase (para AC)
- Gráfico dinâmico da corrente \( I(t) \)
- Interface gráfica intuitiva via Matplotlib

## Modelo Matemático

A equação diferencial do circuito RLC série é dada por:

\[
L \frac{d^2q}{dt^2} + R \frac{dq}{dt} + \frac{q}{C} = V(t)
\]

Onde:
- \( q(t) \): carga no capacitor
- \( I(t) = \frac{dq}{dt} \): corrente

Transformando em um sistema de equações de primeira ordem:

\[
\begin{cases}
\frac{dq}{dt} = I \\
\frac{dI}{dt} = \frac{1}{L} \left( V(t) - R I - \frac{q}{C} \right)
\end{cases}
\]

## Como executar

1. Clone o repositório:
   ```bash
   git clone https://github.com/seu-usuario/simulador-rlc-rungekutta.git
   cd simulador-rlc-rungekutta

2. Instale as dependências:
pip install -r requirements.txt

3. Execute o simulador:
python rlc_simulador.py


Requisitos:
Python 3.x
- numpy
- matplotlib
- ipywidgets (opcional, se usar Jupyter)