# Modelo SIR para Influenza

Author: Querem Hapuque Moura de Lima
Institution: Federal Rural University of Pernambuco, Department of Statistics and Informatics, Brazil
Email: querem.lima@ufrpe.br
#

Uma implementação simples e funcional do modelo compartimental SIR (Susceptible-Infected-Recovered) para análise da dinâmica de transmissão da influenza.

## Descrição

Este projeto implementa o modelo epidemiológico SIR clássico especificamente parametrizado para influenza. O modelo divide a população em três compartimentos e utiliza equações diferenciais para simular a evolução temporal de uma epidemia.

## Requisitos

- Python 3.6+
- NumPy
- Matplotlib
- SciPy

## Instalação

```bash
pip install numpy matplotlib scipy
```

## Uso

Execute o script principal:

```bash
python sir_influenza.py
```

O programa irá:
1. Simular dois cenários epidemiológicos (influenza sazonal vs cepa mais transmissível)
2. Exibir métricas-chave como R₀, tempo de pico, taxa de ataque
3. Gerar gráficos comparativos
4. Realizar análise de sensibilidade para diferentes valores de R₀

## Estrutura do Código

### Classe `SIRModel`

A classe principal que implementa o modelo SIR:

- `__init__(beta, gamma, N)`: Inicializa o modelo com parâmetros epidemiológicos
- `deriv(y, t)`: Define o sistema de equações diferenciais
- `simulate(S0, I0, R0, t)`: Executa a simulação temporal
- `plot_simulation()`: Visualiza os resultados
- `get_peak_info()`: Calcula informações sobre o pico epidêmico
- `calculate_attack_rate()`: Calcula a taxa de ataque final

### Parâmetros do Modelo

- **β (beta)**: Taxa de transmissão (contatos infecciosos por dia)
- **γ (gamma)**: Taxa de recuperação (1/duração da infecção)
- **N**: População total
- **R₀**: Número básico de reprodução (β/γ)

### Cenários Implementados

1. **Influenza Sazonal**: β=0.3, γ=0.1 (R₀=1.5)
2. **Cepa Mais Transmissível**: β=0.5, γ=0.1 (R₀=2.5)

## Saídas

O programa gera:
- Métricas epidemiológicas resumidas no terminal
- Gráficos comparativos da evolução temporal dos compartimentos
- Análise de sensibilidade mostrando relação entre R₀ e taxa de ataque

## Exemplo de Uso Programático

```python
from sir_influenza import SIRModel
import numpy as np

# Criar modelo
model = SIRModel(beta=0.4, gamma=0.1, N=100000)

# Definir condições iniciais
S0, I0, R0 = 99900, 100, 0
t = np.linspace(0, 365, 365)

# Simular
S, I, R = model.simulate(S0, I0, R0, t)

# Analisar resultados
peak_time, peak_infections = model.get_peak_info(I, t)
print(f"Pico: {peak_infections:.0f} casos no dia {peak_time:.0f}")
```