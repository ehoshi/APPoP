CXX=mpicxx
CPPFLAGS=-I/opt/local/include
CXXFLAGS=-Wall -pedantic -Wextra -O2 -Wno-long-long -O0 -g
LDFLAGS=-L/opt/local/lib
LDLIBS= -larmadillo

OBJS=frame.o frm_FXN.o frm_SIM.o frm_MAS.o frm_WRK.o

frame: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

frm_FXN: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

frm_SIM: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

frm_MAS: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

frm_WRK: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)


.PHONY: clean 
clean:
	rm -f frame $(OBJS)

