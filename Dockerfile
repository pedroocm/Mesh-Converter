FROM pymesh/pymesh:latest

WORKDIR /files

ENTRYPOINT ["python", "build.py"]
#ENTRYPOINT ["bash"]
