import os

os.chdir(os.getcwd())
# remove all files end in .lck
for file in os.listdir():
    if file.endswith(".lck"):
        os.remove(file)