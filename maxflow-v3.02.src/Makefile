libmaxflow.so: graph.cpp maxflow.cpp capi.cpp
	g++ -shared -fPIC -Ofast -mtune=native -march=native -o$@ $^

tst: graph.cpp  maxflow.cpp  tst.cpp
	g++ -o$@ $^
