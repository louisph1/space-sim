compile:
	gcc main.c glad/src/glad.c -ospaaaace -O2 -ffast-math -lcglm -lglfw -lGL -lX11 -lpthread -ldl -lm -Iglad/include
	