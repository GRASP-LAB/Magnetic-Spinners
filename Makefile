TARGET = continuous_model

CXX = mpic++  

CXXFLAGS = -std=c++20  -fopenmp -O3 -lm  -fopenmp 

SOURCES = continuous_model_function.cpp continuous_model_physics.cpp continuous_model_main.cpp continuous_model_um.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

HEADERS = continuous_model_function.h continuous_model_physics.h continuous_model_um.h

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)  

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@  

clean:
	rm -f $(OBJECTS) $(TARGET)

run: $(TARGET)
	mpirun -np 4 ./$(TARGET)  