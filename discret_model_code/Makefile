# Compilateur C++
CXX := g++

# Options de compilation
CXXFLAGS := -std=c++20 -fopenmp -O3  

# Options de l'édition des liens
LDFLAGS := -lm

# Nom de l'exécutable
TARGET := discret_model

# Liste des fichiers sources C++
SRCS_CPP := discret_model_main.cpp discret_model_function.cpp
OBJS_CPP := $(SRCS_CPP:.cpp=.o)

# Liste des fichiers d'en-tête
HEADERS := discret_model_function.h

# Règle par défaut (première règle dans le fichier)
all: $(TARGET)

# Règle pour la cible principale
$(TARGET): $(OBJS_CPP)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS_CPP) $(LDFLAGS)

# Règle pour les fichiers sources C++
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règle pour nettoyer les fichiers objets
clean:
	rm -f $(OBJS_CPP) $(TARGET)
