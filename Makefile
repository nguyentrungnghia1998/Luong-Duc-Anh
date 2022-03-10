TARGET= LBM
CXX= g++
CXXFLAGS=-Wall -O3 -fopenmp -lm
SRCDIR= ./src
OBJDIR= ./obj
SOURCE= $(wildcard *.cpp)
OBJECT := $(SOURCE:.cpp=.o)
OUTPUT= out.$(TARGET).txt field*.vtr stop

$(TARGET): $(OBJECT)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECT)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< 

.PHONY: clean
clean: 
	rm -rf $(TARGET) $(OBJECT) $(OUTPUT) $(OBJDIR)