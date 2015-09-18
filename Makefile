metaplot: RunMetaplot.o Metaplot.o
	g++ -Wall -o $@ RunMetaplot.o Metaplot.o
Run.o: RunMetaplot.cpp Metaplot.h
	g++ -Wall -c RunMetaplot.cpp
Metaplot.o: Metaplot.cpp Metaplot.h
	g++ -Wall -c Metaplot.cpp
clean: 
	rm -f Metaplot.o RunMetaplot.o metaplot
