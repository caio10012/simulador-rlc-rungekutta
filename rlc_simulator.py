import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.image as mpimg
import os

# --- Funções da fonte ---
def fonte_dc(t, v):
    return v

def fonte_ac(t, amp, freq):
    return amp * np.sin(2 * np.pi * freq * t)

# --- Equações diferenciais ---
def sistema(t, q, i, R, L, C, fonte, args_fonte):
    V = fonte(t, *args_fonte)
    dqdt = i
    didt = (V - R * i - q / C) / L
    return dqdt, didt

# --- Método de Runge-Kutta 4ª ordem ---
def runge_kutta(t0, q0, i0, dt, t_max, R, L, C, fonte, args_fonte):
    t = [t0]
    q = [q0]
    i = [i0]

    while t[-1] < t_max:
        dq1, di1 = sistema(t[-1], q[-1], i[-1], R, L, C, fonte, args_fonte)
        dq2, di2 = sistema(t[-1]+dt/2, q[-1]+dt*dq1/2, i[-1]+dt*di1/2, R, L, C, fonte, args_fonte)
        dq3, di3 = sistema(t[-1]+dt/2, q[-1]+dt*dq2/2, i[-1]+dt*di2/2, R, L, C, fonte, args_fonte)
        dq4, di4 = sistema(t[-1]+dt, q[-1]+dt*dq3, i[-1]+dt*di3, R, L, C, fonte, args_fonte)

        q.append(q[-1] + dt * (dq1 + 2*dq2 + 2*dq3 + dq4) / 6)
        i.append(i[-1] + dt * (di1 + 2*di2 + 2*di3 + di4) / 6)
        t.append(t[-1] + dt)

    return np.array(t), np.array(i)

# --- Função que atualiza o gráfico ---
def atualizar(val):
    R = slider_R.val
    L = slider_L.val
    C = slider_C.val
    modo = radio_fonte.value_selected

    if modo == "DC":
        V = slider_DC.val
        fonte = fonte_dc
        args = (V,)
    else:
        amp = slider_amp.val
        freq = slider_freq.val
        fonte = fonte_ac
        args = (amp, freq)

    t, i = runge_kutta(0, 0, 0, 0.001, 1.0, R, L, C, fonte, args)
    linha.set_ydata(i)
    linha.set_xdata(t)
    ax_plot.relim()
    ax_plot.autoscale_view()
    fig.canvas.draw_idle()

# --- Interface ---
fig = plt.figure(figsize=(10, 8))
gs = fig.add_gridspec(3, 2)

# Imagem do circuito (se existir)
ax_img = fig.add_subplot(gs[0, 0])
img_path = os.path.join("img", "circuito.png")
if os.path.exists(img_path):
    img = mpimg.imread(img_path)
    ax_img.imshow(img)
    ax_img.axis('off')
    ax_img.set_title("Circuito RLC")
else:
    ax_img.axis('off')
    ax_img.set_title("Imagem não encontrada")

# Gráfico
ax_plot = fig.add_subplot(gs[0, 1])
ax_plot.set_title("Corrente I(t)")
ax_plot.set_xlabel("Tempo (s)")
ax_plot.set_ylabel("Corrente (A)")
ax_plot.grid(True)

# Simulação inicial com tempo de 1.0s
t_ini, i_ini = runge_kutta(0, 0, 0, 0.001, 1.0, 10.0, 0.1, 0.01, fonte_dc, (5.0,))
linha, = ax_plot.plot(t_ini, i_ini, lw=2)

# Sliders
axcolor = 'lightgoldenrodyellow'
ax_R     = plt.axes([0.15, 0.25, 0.65, 0.03], facecolor=axcolor)
ax_L     = plt.axes([0.15, 0.21, 0.65, 0.03], facecolor=axcolor)
ax_C     = plt.axes([0.15, 0.17, 0.65, 0.03], facecolor=axcolor)
ax_DC    = plt.axes([0.15, 0.13, 0.65, 0.03], facecolor=axcolor)
ax_amp   = plt.axes([0.15, 0.09, 0.65, 0.03], facecolor=axcolor)
ax_freq  = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)

slider_R     = Slider(ax_R, 'R (Ohm)', 0.1, 100.0, valinit=10.0)
slider_L     = Slider(ax_L, 'L (H)', 0.001, 10.0, valinit=0.1)
slider_C     = Slider(ax_C, 'C (F)', 0.0001, 1.0, valinit=0.01)
slider_DC    = Slider(ax_DC, 'V (DC)', 0.0, 10.0, valinit=5.0)
slider_amp   = Slider(ax_amp, 'Amplitude (AC)', 0.0, 10.0, valinit=5.0)
slider_freq  = Slider(ax_freq, 'Frequência (Hz)', 1.0, 1000.0, valinit=60.0)

# Tipo de fonte - novo local (topo esquerdo)
ax_radio = plt.axes([0.02, 0.75, 0.1, 0.1], facecolor=axcolor)
radio_fonte = RadioButtons(ax_radio, ('DC', 'AC'), active=0)

# Botão de simular - mais próximo dos sliders
ax_botao = plt.axes([0.83, 0.13, 0.1, 0.04])
botao = Button(ax_botao, 'Simular')
botao.on_clicked(atualizar)

# Inicia com os valores atuais
atualizar(None)

plt.tight_layout()
plt.show()
