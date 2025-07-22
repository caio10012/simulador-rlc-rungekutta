# Importação de bibliotecas essenciais
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.image import imread
from matplotlib.patches import Rectangle
import os
import warnings

# Configurações visuais e de alerta
warnings.filterwarnings('ignore', category=RuntimeWarning)
plt.style.use('ggplot')

# --- Função para carregar imagem do circuito ---
def carregar_imagem_circuito():
    img_path = os.path.join('img', 'circuito.png')
    if os.path.exists(img_path):
        try:
            return imread(img_path)  # Lê a imagem se existir
        except:
            print("Erro ao carregar a imagem. Verifique se o arquivo é uma imagem válida.")
            return None
    else:
        print(f"Arquivo de imagem não encontrado em: {img_path}")
        return None

# Tenta carregar a imagem ao iniciar
circuito_img = carregar_imagem_circuito()

# --- Fontes disponíveis ---
def fonte_dc(t, v):
    # Fonte de tensão contínua (DC)
    return v * np.ones_like(t) if isinstance(t, np.ndarray) else v

def fonte_ac(t, amp, freq):
    # Fonte de tensão alternada senoidal (AC)
    return amp * np.sin(2 * np.pi * freq * t)

# --- Sistema de equações diferenciais ---
def sistema(t, q, i, R, L, C, fonte, args_fonte):
    try:
        V = fonte(t, *args_fonte)  # Tensão da fonte
        dqdt = i  # Derivada da carga é a corrente
        didt = (V - R * i - q / C) / L if L > 1e-10 else 0  # Derivada da corrente
        return dqdt, didt
    except:
        return 0, 0  # Em caso de erro, retorna derivadas nulas

# --- Método de Runge-Kutta de 4ª ordem ---
def runge_kutta(t0, q0, i0, t_max, R, L, C, fonte, args_fonte):
    # Determina passo baseado na frequência da fonte (maior frequência = menor passo)
    if fonte == fonte_ac and args_fonte[1] > 0:
        freq = args_fonte[1]
        dt = 1/(freq * 100)
    else:
        dt = 0.001  # Passo padrão
    
    num_steps = int(t_max / dt) + 1  # Número de pontos
    t = np.linspace(t0, t_max, num_steps)
    q = np.zeros(num_steps)  # Vetor da carga
    i = np.zeros(num_steps)  # Vetor da corrente
    
    q[0] = q0
    i[0] = i0
    
    # Laço principal do método RK4
    for n in range(num_steps - 1):
        current_t = t[n]
        k1_q, k1_i = sistema(current_t, q[n], i[n], R, L, C, fonte, args_fonte)
        k2_q, k2_i = sistema(current_t + dt/2, q[n] + dt*k1_q/2, i[n] + dt*k1_i/2, R, L, C, fonte, args_fonte)
        k3_q, k3_i = sistema(current_t + dt/2, q[n] + dt*k2_q/2, i[n] + dt*k2_i/2, R, L, C, fonte, args_fonte)
        k4_q, k4_i = sistema(current_t + dt, q[n] + dt*k3_q, i[n] + dt*k3_i, R, L, C, fonte, args_fonte)
        q[n+1] = q[n] + dt * (k1_q + 2*k2_q + 2*k3_q + k4_q) / 6
        i[n+1] = i[n] + dt * (k1_i + 2*k2_i + 2*k3_i + k4_i) / 6
    
    return t, i  # Retorna tempo e corrente

# --- Classe para armazenar estado da simulação ---
class CircuitoState:
    def __init__(self):
        self.freeze_updates = False  # Evita múltiplas atualizações simultâneas
        self.last_params = None
        self.update_count = 0  # Contador de atualizações
        self.freq_resonancia = 0
        self.z_impedancia = 0
        self.angulo_fase = 0

# --- Cálculos físicos do circuito ---
def calcular_propriedades(R, L, C, freq):
    # Frequência de ressonância
    if L > 0 and C > 0:
        freq_res = 1/(2 * np.pi * np.sqrt(L * C))
    else:
        freq_res = 0

    # Impedância total e ângulo de fase
    if freq > 0:
        XL = 2 * np.pi * freq * L
        XC = 1/(2 * np.pi * freq * C) if C > 0 else 0
        Z = np.sqrt(R**2 + (XL - XC)**2)
        fase = np.arctan2((XL - XC), R) * 180/np.pi
    else:
        Z = R
        fase = 0
    
    return freq_res, Z, fase

# --- Função que atualiza o gráfico com novos valores ---
def atualizar(val):
    if state.freeze_updates:
        return
    
    try:
        state.freeze_updates = True
        
        # Obtém valores dos sliders
        R = slider_R.val
        L = slider_L.val
        C = slider_C.val
        modo = radio_fonte.value_selected
        
        # Configura a fonte conforme o modo (DC ou AC)
        if modo == "DC":
            V = slider_DC.val
            fonte = fonte_dc
            args = (V,)
            slider_amp.ax.set_visible(False)
            slider_freq.ax.set_visible(False)
            slider_DC.ax.set_visible(True)
            t_max = 1.0
            freq = 0
        else:
            amp = slider_amp.val
            freq = slider_freq.val
            fonte = fonte_ac
            args = (amp, freq)
            slider_amp.ax.set_visible(True)
            slider_freq.ax.set_visible(True)
            slider_DC.ax.set_visible(False)
            t_max = max(1.0, 3/freq) if freq > 0 else 1.0

        # Calcula propriedades do circuito
        state.freq_resonancia, state.z_impedancia, state.angulo_fase = calcular_propriedades(R, L, C, freq)
        
        # Executa simulação
        t, i = runge_kutta(0, 0, 0, t_max, R, L, C, fonte, args)
        
        # Atualiza gráfico
        linha.set_data(t, i)
        max_i = max(np.abs(i)) * 1.2 if len(i) > 0 else 1
        ax_plot.set_ylim(-max_i, max_i)
        ax_plot.set_xlim(0, t_max)
        
        # Atualiza texto de informações
        state.update_count += 1
        update_text.set_text(f"Atualização: {state.update_count}")
        info_text.set_text(
            f"Freq. Ressonância: {state.freq_resonancia:.2f} Hz\n"
            f"Impedância: {state.z_impedancia:.2f} Ω\n"
            f"Ângulo de Fase: {state.angulo_fase:.1f}°\n"
            f"Configuração: R={R:.1f}Ω, L={L:.3f}H, C={C:.4f}F\n"
            f"Fonte: {modo} {'V='+str(args[0])+'V' if modo=='DC' else f'Amp={args[0]}V, Freq={args[1]}Hz'}"
        )
        
        # Realce visual se estiver próximo da ressonância
        if modo == "AC" and abs(freq - state.freq_resonancia) < 0.1 and state.freq_resonancia > 0:
            for patch in ax_plot.patches:
                patch.remove()
            ax_plot.add_patch(Rectangle((0, -max_i), t_max, 2*max_i, 
                             color='yellow', alpha=0.2, label='Ressonância'))
            ax_plot.legend()
        else:
            for patch in ax_plot.patches:
                patch.remove()
        
        fig.canvas.draw_idle()  # Atualiza interface gráfica
        
    finally:
        state.freeze_updates = False

# --- Interface Gráfica com Matplotlib ---
fig = plt.figure(figsize=(16, 10))
fig.subplots_adjust(bottom=0.3, left=0.1, right=0.9, top=0.9)

state = CircuitoState()

# Mostra imagem do circuito (se existir)
if circuito_img is not None:
    ax_img = plt.axes([0.05, 0.7, 0.2, 0.2])
    ax_img.imshow(circuito_img)
    ax_img.axis('off')
    ax_img.set_title("Diagrama do Circuito", pad=10)

# Gráfico principal da corrente
ax_plot = plt.axes([0.3, 0.5, 0.65, 0.4])
ax_plot.set_title("Resposta do Circuito RLC - Corrente I(t)", pad=20)
ax_plot.set_xlabel("Tempo (s)")
ax_plot.set_ylabel("Corrente (A)")
ax_plot.grid(True, alpha=0.3)
ax_plot.axhline(0, color='gray', linestyle='--', alpha=0.5)

# Área com informações da simulação
ax_info = plt.axes([0.3, 0.1, 0.65, 0.3])
ax_info.axis('off')
update_text = ax_info.text(0.02, 0.9, "Atualização: 0", fontsize=10, weight='bold')
info_text = ax_info.text(0.02, 0.5, 
                        "Freq. Ressonância: 0 Hz\n"
                        "Impedância: 0 Ω\n"
                        "Ângulo de Fase: 0°\n"
                        "Configuração: R=0Ω, L=0H, C=0F\n"
                        "Fonte: DC V=0V", 
                        fontsize=10)

# Simulação inicial para preencher gráfico
t_ini, i_ini = runge_kutta(0, 0, 0, 1.0, 10.0, 0.1, 0.01, fonte_dc, (5.0,))
linha, = ax_plot.plot(t_ini, i_ini, lw=2)

# Criação dos sliders
axcolor = 'lightgoldenrodyellow'
slider_x, slider_y, slider_w, slider_h = 0.05, 0.25, 0.2, 0.03
ax_R = plt.axes([slider_x, slider_y, slider_w, slider_h], facecolor=axcolor)
ax_L = plt.axes([slider_x, slider_y-0.05, slider_w, slider_h], facecolor=axcolor)
ax_C = plt.axes([slider_x, slider_y-0.10, slider_w, slider_h], facecolor=axcolor)
ax_DC = plt.axes([slider_x, slider_y-0.15, slider_w, slider_h], facecolor=axcolor)
ax_amp = plt.axes([slider_x, slider_y-0.15, slider_w, slider_h], facecolor=axcolor)
ax_freq = plt.axes([slider_x, slider_y-0.20, slider_w, slider_h], facecolor=axcolor)

slider_R = Slider(ax_R, 'R (Ω)', 0.1, 100.0, valinit=10.0, valstep=0.1)
slider_L = Slider(ax_L, 'L (H)', 0.001, 10.0, valinit=0.1, valstep=0.01)
slider_C = Slider(ax_C, 'C (F)', 0.0001, 1.0, valinit=0.01, valstep=0.001)
slider_DC = Slider(ax_DC, 'V (DC)', 0.0, 20.0, valinit=5.0, valstep=0.1)
slider_amp = Slider(ax_amp, 'Amp (V)', 0.0, 20.0, valinit=5.0, valstep=0.1)
slider_freq = Slider(ax_freq, 'Freq (Hz)', 0.1, 1000.0, valinit=1.0, valstep=0.1)

# Botões para escolher entre fonte DC e AC
ax_radio = plt.axes([0.05, 0.45, 0.2, 0.1], facecolor=axcolor)
radio_fonte = RadioButtons(ax_radio, ('DC', 'AC'), active=0)

# Esconde sliders AC no início
slider_amp.ax.set_visible(False)
slider_freq.ax.set_visible(False)

# Conecta eventos dos sliders à função de atualização
def safe_connect(slider):
    def wrapper(val):
        if not state.freeze_updates:
            atualizar(val)
    slider.on_changed(wrapper)

# Conecta todos os sliders
for slider in [slider_R, slider_L, slider_C, slider_DC, slider_amp, slider_freq]:
    safe_connect(slider)

# Conecta o botão de escolha de fonte
radio_fonte.on_clicked(lambda label: atualizar(None))

# Primeira atualização para carregar tudo
atualizar(None)

# Mostra a interface gráfica interativa
plt.show()
