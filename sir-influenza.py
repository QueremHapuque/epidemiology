import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class SIRModel:
    """
    Modelo compartimental SIR (Susceptible-Infected-Recovered) para influenza
    """
    
    def __init__(self, beta, gamma, N):
        """
        Inicializa o modelo SIR
        
        Parâmetros:
        beta (float): Taxa de transmissão
        gamma (float): Taxa de recuperação
        N (int): População total
        """
        self.beta = beta
        self.gamma = gamma
        self.N = N
        self.R0 = beta / gamma  # Número básico de reprodução
    
    def deriv(self, y, t):
        """
        Sistema de equações diferenciais do modelo SIR
        
        dy/dt = [dS/dt, dI/dt, dR/dt]
        """
        S, I, R = y
        dSdt = -self.beta * S * I / self.N
        dIdt = self.beta * S * I / self.N - self.gamma * I
        dRdt = self.gamma * I
        return dSdt, dIdt, dRdt
    
    def simulate(self, S0, I0, R0, t):
        """
        Simula a evolução temporal do modelo SIR
        
        Parâmetros:
        S0 (float): População inicial susceptível
        I0 (float): População inicial infectada
        R0 (float): População inicial recuperada
        t (array): Array de tempo
        
        Retorna:
        tuple: (S, I, R) arrays com a evolução temporal de cada compartimento
        """
        y0 = [S0, I0, R0]
        sol = odeint(self.deriv, y0, t)
        S, I, R = sol.T
        return S, I, R
    
    def plot_simulation(self, S, I, R, t, title="Evolução da Epidemia de Influenza"):
        """
        Plota os resultados da simulação
        """
        plt.figure(figsize=(10, 6))
        plt.plot(t, S, 'b-', label='Susceptíveis', linewidth=2)
        plt.plot(t, I, 'r-', label='Infectados', linewidth=2)
        plt.plot(t, R, 'g-', label='Recuperados', linewidth=2)
        
        plt.xlabel('Tempo (dias)')
        plt.ylabel('População')
        plt.title(title)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def get_peak_info(self, I, t):
        """
        Calcula informações sobre o pico da epidemia
        """
        peak_idx = np.argmax(I)
        peak_time = t[peak_idx]
        peak_infections = I[peak_idx]
        return peak_time, peak_infections
    
    def calculate_attack_rate(self, S0, S_final):
        """
        Calcula a taxa de ataque (proporção da população que foi infectada)
        """
        return (S0 - S_final) / S0

def main():
    """
    Função principal que executa simulações com diferentes cenários
    """
    
    # Parâmetros baseados na influenza sazonal
    N = 100000  # População total
    
    # Cenário 1: Influenza sazonal típica
    beta1 = 0.3    # Taxa de transmissão (contatos infecciosos por dia)
    gamma1 = 0.1   # Taxa de recuperação (1/duração_infecção)
    
    # Cenário 2: Cepa mais transmissível
    beta2 = 0.5
    gamma2 = 0.1
    
    # Condições iniciais
    I0 = 100      # Casos iniciais
    R0 = 0        # Nenhum recuperado inicialmente
    S0 = N - I0 - R0  # Resto da população é susceptível
    
    # Tempo de simulação (1 ano)
    t = np.linspace(0, 365, 365)
    
    print("=== Simulação do Modelo SIR para Influenza ===\n")
    
    # Simulação Cenário 1
    model1 = SIRModel(beta1, gamma1, N)
    S1, I1, R1 = model1.simulate(S0, I0, R0, t)
    
    peak_time1, peak_infections1 = model1.get_peak_info(I1, t)
    attack_rate1 = model1.calculate_attack_rate(S0, S1[-1])
    
    print(f"Cenário 1 - Influenza Sazonal:")
    print(f"R₀ = {model1.R0:.2f}")
    print(f"Pico da epidemia: {peak_infections1:.0f} infectados no dia {peak_time1:.0f}")
    print(f"Taxa de ataque: {attack_rate1:.1%}")
    print(f"Total de casos: {R1[-1]:.0f}\n")
    
    # Simulação Cenário 2
    model2 = SIRModel(beta2, gamma2, N)
    S2, I2, R2 = model2.simulate(S0, I0, R0, t)
    
    peak_time2, peak_infections2 = model2.get_peak_info(I2, t)
    attack_rate2 = model2.calculate_attack_rate(S0, S2[-1])
    
    print(f"Cenário 2 - Cepa Mais Transmissível:")
    print(f"R₀ = {model2.R0:.2f}")
    print(f"Pico da epidemia: {peak_infections2:.0f} infectados no dia {peak_time2:.0f}")
    print(f"Taxa de ataque: {attack_rate2:.1%}")
    print(f"Total de casos: {R2[-1]:.0f}\n")
    
    # Plotar comparação
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 1, 1)
    plt.plot(t, S1, 'b-', label='Susceptíveis', linewidth=2)
    plt.plot(t, I1, 'r-', label='Infectados', linewidth=2)
    plt.plot(t, R1, 'g-', label='Recuperados', linewidth=2)
    plt.title(f'Cenário 1: Influenza Sazonal (R₀ = {model1.R0:.1f})')
    plt.ylabel('População')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 1, 2)
    plt.plot(t, S2, 'b-', label='Susceptíveis', linewidth=2)
    plt.plot(t, I2, 'r-', label='Infectados', linewidth=2)
    plt.plot(t, R2, 'g-', label='Recuperados', linewidth=2)
    plt.title(f'Cenário 2: Cepa Mais Transmissível (R₀ = {model2.R0:.1f})')
    plt.xlabel('Tempo (dias)')
    plt.ylabel('População')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Análise de sensibilidade para R₀
    print("=== Análise de Sensibilidade ===")
    R0_values = np.linspace(1.1, 3.0, 10)
    final_attack_rates = []
    
    for R0_val in R0_values:
        beta_temp = R0_val * gamma1
        model_temp = SIRModel(beta_temp, gamma1, N)
        S_temp, I_temp, R_temp = model_temp.simulate(S0, I0, R0, t)
        attack_rate_temp = model_temp.calculate_attack_rate(S0, S_temp[-1])
        final_attack_rates.append(attack_rate_temp)
    
    plt.figure(figsize=(8, 6))
    plt.plot(R0_values, final_attack_rates, 'o-', linewidth=2, markersize=6)
    plt.xlabel('R₀')
    plt.ylabel('Taxa de Ataque Final')
    plt.title('Relação entre R₀ e Taxa de Ataque')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()