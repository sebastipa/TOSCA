#bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=1 --ntasks-per-node=40
#SBATCH --exclusive                 # no one enters the node
#SBATCH --time=6:00:00             # time limits 24 hours
#SBATCH --error solutionErr         # std-error file
#SBATCH --output solutionLog        # std-output file
#SBATCH --account=rrg-jbrinker      # account number
#SBATCH --mail-type=ALL
#SBATCH --mail-user=colechr@mail.ubc.ca

mpirun ./tosca > toscaLog
