#include"Metaplot.h"
using namespace std;

int main(int argc, char * argv[])
{
	try{
		UserOpts input(argc, argv);
	
		cout << "Begin splitting bed files ..." << endl;
		vector<Bed *> * bedByChrs = input.splitBedFiles();
		cout << "\tdone." << endl;
		cout << "Begin splitting wig file ..." << endl;
		vector<Wig *> wigByChrs = input.splitWigFile();
		cout << "\tdone." << endl;

	// how to efficiently search for chrs? It looks O(n*n) ...
	// do one awful preprocessing step then use only common chrs
	
	/* the difficulty is that stacks keep moving around and getting reversed */

		vector<string> commonChrs = input.commonChroms();

	// now erase non-common chrs? this all should really be done in a preprocessing step that runs the whole file
	// or just skip if not in common_chrs
	// how to compare string sort values?
// end of preprocessing
		MetaplotRegion * region = new MetaplotRegion[input.getBedNumber()];

		for (int i = 0; i < input.getBedNumber(); i++)
			region[i] = MetaplotRegion(input.getMaxWindow());
		
		calculateMetaplot(input, bedByChrs, wigByChrs, commonChrs, region);
	} // try end paren	
	catch(int e)
	{
		if (e == 1)
		{
			cerr << "Exiting!" << endl;
			return 1;
		}
	}
	return 0;
}

