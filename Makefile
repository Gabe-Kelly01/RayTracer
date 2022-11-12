OBJS	= Tuple.o Rgb.o easyppm.o Sphere.o Plane.o Intersection.o Lighting.o RayTracer.o
SOURCE	= Tuple.cpp Rgb.cpp easyppm.c  Sphere.cpp Plane.cpp Intersection.cpp Lighting.cpp RayTracer.cpp
HEADER	= Tuple.h Rgb.h easyppm.h Sphere.h Plane.h Intersection.h Lighting.h
OUT	= tracer
CC	 = g++
FLAGS	 = -g -c -Wall

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

run: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)
	./$(OUT)

Tuple.o: Tuple.cpp
	$(CC) $(FLAGS) Tuple.cpp -std=c++17

Rgb.o: Rgb.cpp
	$(CC) $(FLAGS) Rgb.cpp -std=c++17

easyppm.o: easyppm.c
	$(CC) $(FLAGS) easyppm.c -std=c++17

Sphere.o: Sphere.cpp
	$(CC) $(FLAGS) Sphere.cpp -std=c++17

Plane.o: Plane.cpp
	$(CC) $(FLAGS) Plane.cpp -std=c++17

Intersection.o: Intersection.cpp
	$(CC) $(FLAGS) Intersection.cpp -std=c++17

Lighting.o: Lighting.cpp
	$(CC) $(FLAGS) Lighting.cpp -std=c++17

RayTracer.o: RayTracer.cpp
	$(CC) $(FLAGS) RayTracer.cpp -std=c++17


clean:
	rm -f $(OBJS) $(OUT)
	rm image.ppm