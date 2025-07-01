# simulador-rlc-rungekutta
Simulador de circuito RLC em série com integração numérica (Runge-Kutta de 4ª ordem) e interface interativa


# Simulador RLC com Método de Runge-Kutta

Este projeto é um simulador interativo de um circuito RLC série, resolvido numericamente utilizando o **método de Runge–Kutta de 4ª ordem**. A simulação permite visualizar a corrente \( I(t) \) em função do tempo para diferentes configurações do circuito.

## Funcionalidades

- Escolha dos parâmetros do circuito: resistência \( R \), indutância \( L \), capacitância \( C \)
- Escolha entre fonte de corrente contínua (DC) ou alternada (AC senoidal)
- Ajuste da tensão DC, ou da amplitude e frequência da fonte AC
- Gráfico dinâmico da corrente \( I(t) \)
- Método de Runge-Kutta de quarta ordem para solução numérica

## Modelo matemático

A equação diferencial do circuito RLC série é:

\[
L \frac{d^2q}{dt^2} + R \frac{dq}{dt} + \frac{q}{C} = V(t)
\]

Definindo:
- \( q \): carga no capacitor
- \( I = \frac{dq}{dt} \): corrente no circuito

Transformamos em um sistema de equações de primeira ordem:

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