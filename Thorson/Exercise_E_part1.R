

# Dimensions
A = 10
T = 40

# Settings
sigma1 = 0.6
sigma2 = 0.1
mu = 1e6
m = 0.4
theta1 = 3
theta2 = 1

# Globals
ages = seq_len(A)

# Randoms
epsilon_at = matrix(
  rnorm(n = A*T),
  nrow = A,
  ncol = T
)

# Initial year
n_at = matrix( nrow=A, ncol = T )
n_at[,1] = mu * exp( sigma1 * epsilon_at[,1] - m * ages )

# Loop through years
for( t in 2:T ){
  n_at[1,t] = mu * exp( sigma1 * epsilon_at[1,t] - m )
  for( a in 2:A ){
    n_at[a,t] = n_at[a-1,t-1] * exp( sigma2 * epsilon_at[a,t] - m )
  }
}

# selectivity at age
s_a = 1 / (1 + exp(-theta2 * (ages-theta1)))
nprime_at = sweep(
  n_at,
  MARGIN = 1,
  STATS = s_a,
  FUN = "*"
)
