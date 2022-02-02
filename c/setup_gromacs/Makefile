# É}ÉNÉçíËã`
CXX			= g++
CXXFLAGS	= -O3 -lm
LDFLAGS		= 
INCLUDES	= -I.
TARGET		= setup_gromacs
OBJS		= setup_gromacs.o

# ê∂ê¨ãKë•
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f $(TARGET) $(OBJS) .nfs* *~ \#* core

.c.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<


