make: 	src/*.cpp
		g++ -fPIC -c src/*.cpp -I./src/
		g++ -shared -o libalg.so *.o -I./src/
		rm -rf *.o

clean:
		rm -rf *.o *.so