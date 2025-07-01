import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import warnings

# Configurações iniciais
warnings.filterwarnings('ignore', category=RuntimeWarning)
plt.style.use('ggplot')

# --- Funções da fonte ---
def fonte_dc(t, v):
    return v * np.ones_like(t) if isinstance(t, np.ndarray) else v

def fonte_ac(t, amp, freq):
    return amp * np.sin(2 * np.pi * freq * t)

# --- Equações diferenciais ---
def sistema(t, q, i, R, L, C, fonte, args_fonte):
    try:
        V = fonte(t, *args_fonte)
        dqdt = i
        didt = (V - R * i - q / C) / L if L > 1e-10 else 0
        return dqdt, didt
    except:
        return 0, 0

# --- Método de Runge-Kutta 4ª ordem adaptativo ---
def runge_kutta(t0, q0, i0, t_max, R, L, C, fonte, args_fonte):
    # Determina o passo de tempo com base na frequência
    if fonte == fonte_ac and args_fonte[1] > 0:  # Se for AC com frequência > 0
        freq = args_fonte[1]
        dt = 1/(freq * 100)  # 100 pontos por ciclo
    else:
        dt = 0.001  # Passo padrão para DC
    
    num_steps = int(t_max / dt) + 1
    t = np.linspace(t0, t_max, num_steps)
    q = np.zeros(num_steps)
    i = np.zeros(num_steps)
    
    q[0] = q0
    i[0] = i0
    
    for n in range(num_steps - 1):
        current_t = t[n]
        
        # Estágios do RK4
        k1_q, k1_i = sistema(current_t, q[n], i[n], R, L, C, fonte, args_fonte)
        k2_q, k2_i = sistema(current_t + dt/2, q[n] + dt*k1_q/2, i[n] + dt*k1_i/2, R, L, C, fonte, args_fonte)
        k3_q, k3_i = sistema(current_t + dt/2, q[n] + dt*k2_q/2, i[n] + dt*k2_i/2, R, L, C, fonte, args_fonte)
        k4_q, k4_i = sistema(current_t + dt, q[n] + dt*k3_q, i[n] + dt*k3_i, R, L, C, fonte, args_fonte)
        
        q[n+1] = q[n] + dt * (k1_q + 2*k2_q + 2*k3_q + k4_q) / 6
        i[n+1] = i[n] + dt * (k1_i + 2*k2_i + 2*k3_i + k4_i) / 6
    
    return t, i

# --- Classe para gerenciar o estado ---
class CircuitoState:
    def __init__(self):
        self.freeze_updates = False
        self.last_params = None

# --- Função que atualiza o gráfico ---
def atualizar(val):
    if state.freeze_updates:
        return
        
    try:
        state.freeze_updates = True
        
        R = slider_R.val
        L = slider_L.val
        C = slider_C.val
        modo = radio_fonte.value_selected

        if modo == "DC":
            V = slider_DC.val
            fonte = fonte_dc
            args = (V,)
            slider_amp.ax.set_visible(False)
            slider_freq.ax.set_visible(False)
            slider_DC.ax.set_visible(True)
            t_max = 1.0  # 1 segundo para DC
        else:
            amp = slider_amp.val
            freq = slider_freq.val
            fonte = fonte_ac
            args = (amp, freq)
            slider_amp.ax.set_visible(True)
            slider_freq.ax.set_visible(True)
            slider_DC.ax.set_visible(False)
            t_max = max(1.0, 3/freq) if freq > 0 else 1.0  # Mostra pelo menos 3 ciclos

        # Verifica se os parâmetros realmente mudaram
        current_params = (R, L, C, modo, *args)
        if state.last_params == current_params:
            return
        state.last_params = current_params

        t, i = runge_kutta(0, 0, 0, t_max, R, L, C, fonte, args)
        
        # Atualização dos dados
        linha.set_data(t, i)
        
        # Ajuste dinâmico dos eixos
        max_i = max(np.abs(i)) * 1.2 if len(i) > 0 else 1
        ax_plot.set_ylim(-max_i, max_i)
        ax_plot.set_xlim(0, t_max)
        
        fig.canvas.draw_idle()
        
    finally:
        state.freeze_updates = False

# --- Configuração da interface ---
fig = plt.figure(figsize=(12, 8))
fig.subplots_adjust(bottom=0.3, left=0.15)

# Estado global
state = CircuitoState()

# Área do gráfico principal
ax_plot = plt.axes([0.1, 0.4, 0.85, 0.55])
ax_plot.set_title("Resposta do Circuito RLC - Corrente I(t)", pad=20)
ax_plot.set_xlabel("Tempo (s)")
ax_plot.set_ylabel("Corrente (A)")
ax_plot.grid(True, alpha=0.3)
ax_plot.axhline(0, color='gray', linestyle='--', alpha=0.5)

# Simulação inicial
t_ini, i_ini = runge_kutta(0, 0, 0, 1.0, 10.0, 0.1, 0.01, fonte_dc, (5.0,))
linha, = ax_plot.plot(t_ini, i_ini, lw=2)

# Sliders
axcolor = 'lightgoldenrodyellow'
slider_x, slider_y, slider_w, slider_h = 0.25, 0.25, 0.6, 0.03
ax_R = plt.axes([slider_x, slider_y, slider_w, slider_h], facecolor=axcolor)
ax_L = plt.axes([slider_x, slider_y-0.05, slider_w, slider_h], facecolor=axcolor)
ax_C = plt.axes([slider_x, slider_y-0.10, slider_w, slider_h], facecolor=axcolor)
ax_DC = plt.axes([slider_x, slider_y-0.15, slider_w, slider_h], facecolor=axcolor)
ax_amp = plt.axes([slider_x, slider_y-0.15, slider_w, slider_h], facecolor=axcolor)
ax_freq = plt.axes([slider_x, slider_y-0.20, slider_w, slider_h], facecolor=axcolor)

slider_R = Slider(ax_R, 'Resistência (Ω)', 0.1, 100.0, valinit=10.0, valstep=0.1)
slider_L = Slider(ax_L, 'Indutância (H)', 0.001, 10.0, valinit=0.1, valstep=0.01)
slider_C = Slider(ax_C, 'Capacitância (F)', 0.0001, 1.0, valinit=0.01, valstep=0.001)
slider_DC = Slider(ax_DC, 'Tensão DC (V)', 0.0, 20.0, valinit=5.0, valstep=0.1)
slider_amp = Slider(ax_amp, 'Amplitude AC (V)', 0.0, 20.0, valinit=5.0, valstep=0.1)
slider_freq = Slider(ax_freq, 'Frequência (Hz)', 0.1, 1000.0, valinit=1.0, valstep=0.1)

# Radio buttons para seleção do tipo de fonte
ax_radio = plt.axes([0.05, 0.4, 0.08, 0.1], facecolor=axcolor)
radio_fonte = RadioButtons(ax_radio, ('DC', 'AC'), active=0)

# Configuração inicial - mostra DC, esconde AC
slider_amp.ax.set_visible(False)
slider_freq.ax.set_visible(False)

# Conecta os eventos
def safe_connect(slider):
    def wrapper(val):
        if not state.freeze_updates:
            atualizar(val)
    slider.on_changed(wrapper)

for slider in [slider_R, slider_L, slider_C, slider_DC, slider_amp, slider_freq]:
    safe_connect(slider)

radio_fonte.on_clicked(lambda label: atualizar(None))

plt.show()