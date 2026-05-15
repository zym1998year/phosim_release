all:
	 rm -f source/raytrace/*.o; rm -f source/trim/*.o; rm -f source/e2adc/*.o; rm -f source/instrument/*.o; rm -f source/atmosphere/*.o ; rm -f source/ancillary/*.o; rm -f source/catalogtools/*.o
	 bash -c "rm -f source/raytrace/ext.h"
	 bash -c "grep -q EXT bin/setup && echo '#define EXT' >> source/raytrace/ext.h || echo '#undef EXT' >> source/raytrace/ext.h"
	 bash -c ". bin/setup ; cd source/raytrace ; make; mv raytrace ../../bin/."
	 bash -c ". bin/setup ; cd source/trim ; make; mv trim ../../bin/."
	 bash -c ". bin/setup ; cd source/e2adc ; make; mv e2adc ../../bin/."
	 bash -c ". bin/setup ; cd source/instrument ; make; mv instrument ../../bin/."
	 bash -c ". bin/setup ; cd source/atmosphere ; make; mv atmosphere ../../bin/."
	 bash -c ". bin/setup ; cd source/ancillary ; make "
	 bash -c ". bin/setup ; cd source/catalogtools ; make;  mv phosimcatgen  ../../bin/."
	 bash -c ". bin/setup ; cd validation ; make "
	 bash -c "rm -f bin/version"
	 bash -c "if [ -d .git ]; then git log -1 --format=commit' '%H >> bin/version ; fi"
	 bash -c "if [ -d .git ]; then git describe --tags >> bin/version ; fi"
	 bash -c "if [ ! -d .git ]; then echo 'commit none' >> bin/version ; fi"
	 bash -c "if [ ! -d .git ]; then pwd >> bin/version ; fi"
	 @echo ' '
	 @echo 'The Photon Simulator has been compiled successfully!'
	 @echo ' '
	 @echo 'Type "phosim_gui" or "phosim" to begin!'
	 @echo ' '
