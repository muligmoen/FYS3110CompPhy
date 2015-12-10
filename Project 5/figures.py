import subprocess

def get_stdout(cmd):
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
  return p.stdout.read()


cmd = ['./project5', str(dt), str(dx), str(T)]