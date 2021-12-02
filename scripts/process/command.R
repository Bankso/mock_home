#!/opt/conda/envs/SCAR_env/bin/Rscript

# Make functions available

library('SCAR')
library('stringr')

# Pull input variables from Singularity
env_vars <- Sys.getenv(c("func", "args"), names=TRUE)

# Make sure env_vars can be accessed/assigned
func <- env_vars[['func']]
args <- env_vars[['args']]

args <- str_c(args, sep = ' ')
system2(func, args=args, stderr = str_c(func, '_log.txt'))