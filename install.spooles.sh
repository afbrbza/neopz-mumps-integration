sudo mkdir -p /opt/spooles/
cd /opt/spooles/
sudo wget https://www.netlib.org/linalg/spooles/spooles.2.2.tgz
sudo tar -xvf spooles.2.2.tgz
sudo rm -rf spooles.2.2.tgz 

# 1. Editar Make.inc para configurar compilador GCC e otimização O3
sudo sed -i 's/^  CC = \/usr\/lang-4.0\/bin\/cc/  CC = gcc/' Make.inc
sudo sed -i 's/^# CC = gcc/# CC = gcc/' Make.inc
sudo sed -i 's/^  OPTLEVEL = -O$/# OPTLEVEL = -O/' Make.inc
sudo sed -i 's/^# OPTLEVEL = -O3$/  OPTLEVEL = -O3/' Make.inc

# 2. Limpar qualquer build anterior
sudo make clean

# 3. Compilar a biblioteca (gera spooles.a)
sudo make lib

# Resultado esperado:
# 
# - Biblioteca compilada: spooles.a (~2.7 MB)
# - Suporte a threads POSIX habilitado por padrão
# - Otimização: -O3 (release mode)