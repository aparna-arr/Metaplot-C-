#include"Metaplot.h"
using namespace std;

/* User Opts functions */

UserOpts::UserOpts(int argc, char *argv[])
{
	if (argc < 5)
	{
		printUsage();
		throw 1;
	}

	step = -1;
	window = -1;
	mode = -1;
	preprocessBedOpt = false;
	preprocessWigOpt = false;
	splitWigDir = "";
	splitBedDir = "";
	wigsSplit = false;
	bedsSplit = false;
	isMonteCarlo = false;
	bedsFromDir = false;
	bedDir = "";

	cout << "Handling opts" << endl;
	string args = handleOpts(argc, argv);

	stringstream input(args);

	if (!(input >> maxWindow))
	{
		cout << "Error! Your maxWindow is not a valid number" << endl;
		throw 1;
	}

	cout << "Your maxWindow is:" << endl << maxWindow << endl;
	
	input >> wigFile;
	cout << "Your wigfile is:" << endl << wigFile << endl;

	vector<string> tmpBedsAndNames;

	string bedOrName;

	// if bedsFromDir, names are filenames

	if (!bedsFromDir)
	{
		while((input >> bedOrName))
		{
			tmpBedsAndNames.push_back(bedOrName);
		}
	}
	else
	{
		// get all filesnames in bedsDir
		readBedsFromDir(tmpBedsAndNames, bedDir);
	}

	if ( tmpBedsAndNames.size() % 2 != 0 )
	{
		cout << "Error! Something is wrong with your bedfiles and names. Odd number returned." << endl;

		cout << "tmpBedsAndNames.size() is " << tmpBedsAndNames.size() << endl;
		cout << "Your input is:" << endl;

		for (vector<string>::iterator iter = tmpBedsAndNames.begin(); iter != tmpBedsAndNames.end(); iter++)
		{
			cout << *iter << endl;
		}	
		throw 1;

	}

	numBeds = tmpBedsAndNames.size() / 2;

	string bed;	

	cout << "Your bedfiles are:" << endl;

	for (vector<string>::iterator iter = tmpBedsAndNames.begin(); iter != tmpBedsAndNames.end() - numBeds; iter++)
	{
		cout << *iter << endl;
		bedFiles.push_back(*iter);
	}

	string name;
	cout << "Your names are:" << endl;
	
	for (vector<string>::iterator iter = tmpBedsAndNames.begin() + numBeds; iter != tmpBedsAndNames.end(); iter++)
	{
		cout << *iter << endl;
		names.push(*iter);
	}

	cout << "preprocessing steps ... " << endl;
	if (preprocessBedOpt)
	{
		cout << "Begin preprocess bed" << endl;
		preprocessBed();
	}

	if (preprocessWigOpt)
	{
		cout << "Begin preprocess wig" << endl;
		preprocessWig();

	}

	cout << "done" << endl;
/* Not using currently
	if (!validateBeds())
	{
		cout << "Your bedfiles are not validating!" << endl;
		throw 1;
	}

	if (!validateWig())
	{
		cout << "Your wigfile is not validating!" << endl;
		throw 1;
	}

	if (!validateNames())
	{
		cout << "Your names are not validating!" << endl;
		throw 1;
	}
	*/
}

void UserOpts::preprocessBed()
{
	if (mode == -1)
		return;

	for (vector<string>::iterator iter = bedFiles.begin(); iter != bedFiles.end(); iter++)
	{
		ifstream file((*iter).c_str());
		string line;
		ofstream tmpfile("tmp.bed");
		while(getline(file, line))
		{
			int start, end;
			char strand;
			string chr;
			stringstream linestream(line);
	
			linestream >> chr >> start >> end >> strand;
			
			int more, less;

			if (mode == 0) // aka TSS
			{
				if (strand == '+')
				{
					more = start + (maxWindow / 2);
					less = start - (maxWindow / 2);
				}
				else
				{
					more = end + (maxWindow / 2);
					less = end - (maxWindow / 2);
				}
			} // if	
			else if (mode == 1) // centered
			{
				int middle = (end - start) / 2 + start;
				more = middle + (maxWindow / 2);
				less = middle - (maxWindow / 2);
			} // else if	

			tmpfile << chr << "\t" << less << "\t" << more << "\t" << strand << endl;
			
		} // while
		file.close();
		tmpfile.close();
		(*iter) = (*iter) + ".tmp";
		rename("tmp.bed", (*iter).c_str());
	} // for
}


string UserOpts::handleOpts(int argc, char * argv[])
{
	// opts are:
	// --step X
	// --window X
	// --debug
	// --preprocessBed X
	// --preprocessWig 
	// --readSplitWig <DIR>
	// --readSplitBed <DIR>
	// --readAllBedsInDir <DIR> // instead of providing list of bedfiles
	// --monteCarlo
	string allArgs = "";

	for (int i = 0; i < argc ; i++)
	{
		string tmp(argv[i]);
		allArgs += tmp + " ";
	}


	stringstream args(allArgs);
	string sub;

	args >> sub; // $0
	args >> sub; // $0

	int pos;

	while((pos = sub.find("--")) != (signed) string::npos)
	{
		string opt = sub.substr(pos+2, sub.length());
	
		if (opt == "step" || opt == "window" || opt == "preprocessBed")
		{
			int num;
			if (!(args >> num))
			{
				cout << "Opt " << opt << " requires an INT argument." << endl;
				throw 1;
			}
			
			if (opt == "step")
			{
				step = num;
			}
			else if (opt == "window")
			{
				window = num;
			}
			else 
			{
				mode = num;
				preprocessBedOpt = true;
			}
			
		}
		else if (opt == "readSplitWig" || opt == "readSplitBed" || opt == "readAllBedsInDir")
		{
			string dir;
			if (!(args >> dir))
			{
				cout << "Opt " << opt << " requires a STRING argument." << endl;
				throw 1;
			}

			if (opt == "readSplitBed")
			{
				splitBedDir = dir;
				bedsSplit = true;
			}
			else if (opt == "readSplitWig")
			{
				splitWigDir = dir;
				wigsSplit = true;
			}
			else if (opt == "readAllBedsInDir")
			{
				bedDir = dir;
				bedsFromDir = true;
			}
		}
		else if (opt == "preprocessWig" || opt == "monteCarlo")
		{
			if (opt == "preprocessWig")
				preprocessWigOpt = true;
			else if (opt == "monteCarlo")
				isMonteCarlo = true;
		}
		else
		{
			cout << "Unrecognized opt: " << opt << endl;
			throw 1;
		}

		args >> sub;
	}	
	
	string rem;
	getline(args, rem);

	return sub + " " + rem;
}

void UserOpts::preprocessWig(void)
{
	if (window == -1) // no window passed
		window = 100;

	if (step == -1) // no step passed
		step = 1;
		
	string line;
	ifstream file(wigFile.c_str());


	int currSpan = 0;
	string currChr = "INIT";

	vector<Peak> * wigBlocks = new vector<Peak>;	

	while(getline(file, line))
	{
		stringstream linestream(line);
		stringstream test(line);
		int pos;
		double val;
		if (!(test >> pos))
		{
			string substring;
			linestream >> substring; // variableStep
			linestream >> substring; // chrom=

			string chr = substring.substr(substring.find("=") + 1, substring.length());

			linestream >> substring; // span=

			string span = substring.substr(substring.find("=") + 1, substring.length());
			
			if (chr != currChr)
			{
			//	tmpfile << "variableStep\tchrom=" << chr << "\tspan=" << window << endl;

				if (currChr != "INIT")
				{
					calcSlidingWindow(wigBlocks);
					wigBlocks->clear();
				}

				currChr = chr;
			} // if
			stringstream itoS(span);
			itoS >> currSpan;
		}	
		else
		{
			linestream >> pos >> val;
			Peak tmp;
			tmp.start = pos;
			tmp.end = pos + currSpan;
			tmp.chr = currChr;
			tmp.value = val;
			wigBlocks->push_back(tmp);
		}
	}		
	calcSlidingWindow(wigBlocks);
	wigBlocks->clear();
	file.close();
	wigFile = wigFile + ".tmp";
	rename("tmp.wig", wigFile.c_str());

}	

void UserOpts::calcSlidingWindow(vector<Peak> * wigBlocks)
{
	cout << "\tCalculating sliding window for chr " << (wigBlocks->begin())->chr << endl;
	ofstream tmpfile;
	tmpfile.open("tmp.wig", ios::app);
// span = window much better with same trends than span=1
	tmpfile << "variableStep chrom=" << (wigBlocks->begin())->chr << " span=" << window << endl;

	for (vector<Peak>::iterator iter = wigBlocks->begin(); iter < wigBlocks->end(); iter++)
	{
		vector<Peak>::iterator iterBack = iter;
		vector<Peak>::iterator iterForward = iter;
		
		while (iterBack != wigBlocks->begin() && iterBack->end > iter->start - (window - 1))
			iterBack--;

		while(iterForward != wigBlocks->end() && iterForward->start < iter->end + (window - 1))
			iterForward++;

		double runningAvg = -1;
		vector<Peak>::iterator saveIterBack = iterBack;

		for (int pos = (iter->start - (window - 1)); pos < (iter->end + (window - 1)); pos += step)
		{
			double value = 0;
			if (runningAvg < 0)
			{

				while (iterBack != iterForward + 1 && iterBack->start < pos + window)
				{
					value+=iterBack->value;
					iterBack++;
				}
				iterBack--;
	
				runningAvg = value;
			
			}
			else
			{
				if (saveIterBack->end < pos)
					runningAvg -= saveIterBack->value;

				saveIterBack = iterBack;
				while (iterBack != iterForward + 1 && iterBack->start < pos + window)
				{
					if (iterBack->end >= pos)
					{
						runningAvg+=iterBack->value;
					}

					iterBack++;
				}
				
			}
			
			if (runningAvg > 0)
				tmpfile << pos << "\t" << (double)runningAvg/(double)window << endl;
		}
	
		iter = iterBack-1;
	}
	tmpfile.close();
	cerr << "\t\tdone" << endl;
}

void UserOpts::debugVectorSize(void)
{
	cerr << "bedFiles size is " << bedFiles.size() << endl;
}

void UserOpts::printUsage(void)
{
	cout << "usage: metaplot [opts] <max_window> <wigfile> <bedfile1 bedfile2 ...> <bedname1 bedname2 ...>" << endl;
	cout << "opts are:" << endl;
	cout << "\t--window <INT> : default 100. Window size to smooth wig file" << endl;
	cout << "\t--step <INT> : default 1. Step size to smooth wig file" << endl;
	cout << "\t--debug : turn on debug mode" << endl;
	cout << "\t--preprocessBed <INT> : preprocess bed. Modes are:" << endl;
	cout << "\t\t0 : TSS" << endl;
	cout << "\t\t1 : Centered" << endl;
	cout << "\t--preprocessWig : preprocess wig." << endl;
	cout << "\t--readSplitWig <DIR>" << endl;
	cout << "\t--readSplitBed <DIR>" << endl;
	cout << "\t--monteCarlo : run a signal simulation. All regions to be used MUST be input as beds!" << endl;
	cout << "\t--readAllBedsInDir <DIR> : instead of listing bed files on cli, read from dir. Names = bed file names. Useful with --monteCarlo. Do NOT input <bedfiles> and <bednames> on cli.";
}	

int UserOpts::getMaxWindow(void)
{
	return maxWindow;
}

bool UserOpts::validateBeds(void)
{
	return true;	
}

bool UserOpts::validateWig(void)
{
	return true;
}

bool UserOpts::validateNames(void)
{
	return true;
}

bool UserOpts::monteCarlo(void)
{
	if (isMonteCarlo)
		return true;;
	
	return false;
}

void UserOpts::readBedsFromDir(vector<string> & tmpNamesAndBeds, string dir)
{
	stack<string> files = getFilesInDir(dir);

	vector<string> names;

	while(!files.empty())
	{
		tmpNamesAndBeds.push_back(files.top());
		names.push_back(files.top());
		files.pop();			
	}

	tmpNamesAndBeds.insert(tmpNamesAndBeds.end(), names.begin(), names.end());
}

vector<Bed *> * UserOpts::splitBedFiles(void)
{
	vector<Bed *> * arOfBedfiles = new vector<Bed*>[numBeds];
	int currIndex = 0;

	if (bedsSplit)
	{
		// beds already split
		// get filenames and chrNames and bed number only from files
		// NOTE: assumes previous run's names & beds were in the EXACT SAME ORDER as this time's! VERY IMPORTANT
		return readBedSplit(arOfBedfiles);
	}

//	cerr << "bedFiles len " << bedFiles.end() - bedFiles.begin() << endl;
//	cerr <<"bedFiles size " << bedFiles.size() << endl;


	for (vector<string>::iterator bedIter = bedFiles.begin(); bedIter != bedFiles.end(); bedIter++)
	{
		cout << "Bedfile : " << *bedIter << endl;
		ifstream fstream((*bedIter).c_str());
		string line;
		
		string prevChr = "INIT";
		
		Bed * currBed = new Bed(); 

		while(getline(fstream, line))
		{
			stringstream linestream(line);

			string chrom;
			int start, end;
			char strand;

			linestream >> chrom;
			linestream >> start >> end;
			linestream >> strand;

			if (chrom != prevChr)
			{
				cout << "\tOn chromosome " << chrom << endl;
				if (prevChr != "INIT")
				{
					currBed->printPeaks();
					Bed * tmp = new Bed();
					*tmp = *currBed;
					arOfBedfiles[currIndex].push_back(tmp);
				}
				
				*currBed = Bed(chrom, maxWindow);
				currBed->generateFilename(bedIter - bedFiles.begin());
				prevChr = chrom;

				if (bedIter - bedFiles.begin() == 0)
				{
					// add every chr to the common chrs list
					Freq newChr;	
					newChr.chr = chrom;
					newChr.freq = 0;
					chroms.push_back(newChr);
				}		
				else
				{
					// search + increment
					// ASSUMING we only hit this conditional once per chr
					for (vector<Freq>::iterator iter = chroms.begin(); iter != chroms.end(); iter++)
					{
						if (chrom == (*iter).chr)
						{
							(*iter).freq++;
							break;
						}
					}
				}
			}
			currBed->addPeak(start, end, strand);
		}

		currBed->printPeaks();
		Bed * tmp = new Bed();
		*tmp = *currBed;
		arOfBedfiles[currIndex].push_back(tmp);

		fstream.close();
		currIndex++;
	}	
	return arOfBedfiles;
}

vector<Wig *> UserOpts::splitWigFile(void)
{
	ifstream fstream(wigFile.c_str());
	string line;

	Wig * currWig = new Wig();
	int currSpan = 0;
	string currChr = "INIT";
	vector<Wig *> wigSplit;

	if (wigsSplit)
	{
		// wig already split
		// get filenames and chrNames only from files
		return readWigSplit(wigSplit);
	}

	while(getline(fstream, line))
	{
		stringstream linestream(line);
		stringstream test(line);
		int pos;
		double val;
		string substring;
		if (!(test >> pos))
		{
			// on a non-integer line
			// FIXME do a check here for variableStep, track lines, and #'s
			linestream >> substring; // variableStep
			linestream >> substring; // chrom=
	
			string chr = substring.substr(substring.find("=") + 1, substring.length());
			
			linestream >> substring; // span=

			string span = substring.substr(substring.find("=") + 1, substring.length());

			if (chr != currChr)
			{
				cout << "\tOn chromosome " << chr << endl;
				if (currChr != "INIT")
				{
					// print prev Chr peaks
					currWig->printPeaks();
					Wig * tmp = new Wig();
					*tmp = *currWig;
					wigSplit.push_back(tmp);	
				}
				
				currChr = chr;
				*currWig = Wig(chr);
				currWig->generateFilename();
				
				for (vector<Freq>::iterator iter = chroms.begin(); iter != chroms.end(); iter++)
				{
					if ((*iter).chr == chr)
					{
						(*iter).freq++;
						break;
					}
				}
			}
			stringstream itoS(span);
			itoS >> currSpan;
		}
		else
		{
			linestream >> pos >> val;
			currWig->addPeak(pos, pos+currSpan, val);
		}
	}

	currWig->printPeaks();
	Wig * tmp = new Wig();
	*tmp = *currWig;
	wigSplit.push_back(tmp);	

	fstream.close();

	return wigSplit;
}

vector<Bed *> * UserOpts::readBedSplit(vector<Bed *> * &arOfBedfiles)
{
	stack<string> splitFiles = getFilesInDir(splitBedDir);
	int pos;
	while(!splitFiles.empty())
	{
		if ( (pos = (splitFiles.top()).rfind("_")) != (signed) string::npos )
		{
			string chr = (splitFiles.top()).substr(0,pos);
			chr = chr.substr(chr.rfind("/")+1);
			stringstream findBedNum((splitFiles.top()).substr(pos+1));
			int bedNum;

			if (!(findBedNum >> bedNum))
			{
				cout << "There is no bednum after _ in file " << splitFiles.top() << endl;
				throw 1;

			}

			Bed * tmp = new Bed(chr, maxWindow);
			tmp->setFilename(splitFiles.top());
			
			// this is prob not as efficient as SplitBedFiles 
			// because this does not intersect chrs in bed files, but merges them
			// because no other way
			// so vector is bigger but other than that should be okay

			if (chroms.empty())
			{
				Freq newChr;	
				newChr.chr = chr;
				newChr.freq = 0;
				chroms.push_back(newChr);
			}	
			else
			{
				bool found = false;
				for (vector<Freq>::iterator iter = chroms.begin(); iter != chroms.end(); iter++)
				{
					if (chr == (*iter).chr)
					{
						iter->freq++;
						found = true;
						break;
					}
				}

				if (!found)
				{
					Freq newChr;	
					newChr.chr = chr;
					newChr.freq = 0;
					chroms.push_back(newChr);
				}
			}

			arOfBedfiles[bedNum].push_back(tmp);	

		}
		else
		{
			cout << "Split file " << splitFiles.top() << " has a weird filename (not *_*)!" << endl;	
			throw 1;

		}
		splitFiles.pop();
	}

	return arOfBedfiles;
}

vector<Wig *> UserOpts::readWigSplit(vector<Wig *> &wigSplit)
{
	stack<string> splitFiles = getFilesInDir(splitWigDir);
	int pos;
	while(!splitFiles.empty())
	{
		if ( (pos = (splitFiles.top()).find("_wig.bed")) != (signed) string::npos )
		{
			string chr = (splitFiles.top()).substr(0,pos);
			chr = chr.substr(chr.rfind("/")+1);
			Wig * tmp = new Wig(chr);
			tmp->setFilename(splitFiles.top());
		

		// NOTE: ALWAYS run splitBedFiles() before splitWigFiles() or else all chroms get messed up 	
			for (vector<Freq>::iterator iter = chroms.begin(); iter != chroms.end(); iter++)
			{
				if (iter->chr == chr)
				{
					iter->freq++;			
					break;
				}
			}
	
			wigSplit.push_back(tmp);

		}
		else
		{
			cout << "Split file " << splitFiles.top() << " has a weird filename (not *_wig.bed)!" << endl;	
			throw 1;
		}
		
		splitFiles.pop();
	}

	return wigSplit;
}

stack<string> UserOpts::getFilesInDir(string dir)
{
	stack<string> filenames;
	
	DIR *dp;
	struct dirent *dirp;
	struct stat filestat;
	string filepath;

	dp = opendir(dir.c_str());
	
	if (dp == NULL)
	{
		cout << "Error opening dir " << dir << endl;
		throw 1;
	}

	while ((dirp = readdir(dp)))
	{
		filepath = dir + "/" + dirp->d_name;
		// skip dirs or weird files
		if (stat(filepath.c_str(), &filestat)) continue;
		if (S_ISDIR(filestat.st_mode)) continue;
		filenames.push(filepath);
	}

	return filenames;
}

int UserOpts::getBedNumber(void)
{
	return numBeds;
}

vector<string> UserOpts::commonChroms(void)
{
	vector<string> commonChrs;
	for (vector<Freq>::iterator iter = chroms.begin(); iter != chroms.end(); iter++) 
	{
		if ((*iter).freq == numBeds)
			commonChrs.push_back((*iter).chr);
	}

	return commonChrs;
}

string UserOpts::getNameString(void)
{
//	cerr << "NAMESTR called" << endl;	
	string namestr = "";

//	if (names.empty())
//		cerr << "names is empty!" << endl;

	while(!names.empty())
	{
		namestr = names.top() + "\t" + namestr;
//		cerr << "adding " << names.top() << " to namesR" << endl;
		namesR.push(names.top());
		names.pop();
	}

	return namestr;
}	

string UserOpts::getNameStringR(void)
{
//	cerr << "beginning of getNameStringR" << endl;

//	if (namesR.empty())
//		cerr << "namesR is empty" << endl;

	string namestr = "'" + namesR.top() + "'";

//	cerr << "\tbefore pop()" << endl;

	namesR.pop();
	
//	cerr << "\tbefore while" << endl;
	
	while(!namesR.empty())
	{
		namestr = "'" + namesR.top() + "', " + namestr;
		namesR.pop();
	}
//	cerr << "end of getNameStringR" << endl;

	return namestr;
	
}
/* End User Opts functions */

/* Chromosome functions */


void Chromosome::clearPeaks(void)
{
	while(!peaks.empty())
	{
		peaks.pop();
	}
}

void Chromosome::setFilename(string filename)
{
	fileName = filename;
}

string Chromosome::getChr(void)
{
	return chrName;
}

void Chromosome::addPeak(void) 
{

}

void Chromosome::printPeaks(void)
{

}

void Chromosome::readPeaks(void)
{

}

Peak * Chromosome::getCurrPeak(void)
{
	if (peaks.empty())
		return NULL;

	return &(peaks.top());
}

void Chromosome::nextPeak(void)
{
	if (!peaks.empty())
		peaks.pop();
}

void Chromosome::generateFilename(void)
{

}

/* End Chromosome functions */

/* Wig functions */

void Wig::unstack(void)
{
	unstacked.push(peaks.top());
	peaks.pop();
}

void Wig::restack(void)
{
	while(!unstacked.empty())
	{
		peaks.push(unstacked.top());
		unstacked.pop();
	}
}

void Wig::printPeaks(void)
{
	ofstream fstream(fileName.c_str());
	while(!peaks.empty())
	{
		Peak tmp = peaks.top();
		fstream << chrName << "\t" << tmp.start << "\t" << tmp.end << "\t" << tmp.value << endl;
		peaks.pop();
	}
	fstream.close();		
}

void Wig::readPeaks(void)
{
	ifstream fstream(fileName.c_str());
	string line;

	while(getline(fstream, line))
	{
		stringstream linestream(line);
		int start, end;
		double value;
		string chrom;
		
		linestream >> chrom >> start >> end >> value;

		addPeak(start, end, value);
	}	
	
	fstream.close();

}

void Wig::addPeak(int start, int end, double value)
{
	Peak tmp;

	tmp.start = start;
	tmp.end = end;
	tmp.value = value;

	peaks.push(tmp);

}

void Wig::generateFilename(void)
{
	fileName = chrName + "_wig" + ".bed";
}

/* End Wig functions */

/* Bed functions */

void Bed::printPeaks(void)
{
	ofstream fstream(fileName.c_str());
	while(!peaks.empty())
	{
		Peak tmp = peaks.top();
		fstream << chrName << "\t" << tmp.start << "\t" << tmp.end << "\t" << tmp.strand << endl;
		peaks.pop();
	}
	fstream.close();		
}

void Bed::readPeaks(void)
{
	ifstream fstream(fileName.c_str());
	string line;

	while(getline(fstream, line))
	{
		stringstream linestream(line);
		int start, end;
		char strand;
		string chrom;
		
		linestream >> chrom >> start >> end >> strand;

		addPeak(start, end, strand);
	}	
	
	fstream.close();
}

void Bed::addPeak(int start, int end, char strand)
{
	Peak tmp;

	tmp.start = start;
	tmp.end = end;
	tmp.strand = strand;

	peaks.push(tmp);
}

void Bed::generateFilename(int filenum)
{
	stringstream intToString;
	intToString << filenum;
	fileName = chrName + "_" + intToString.str() + ".bed";
}

/* End Bed functions */

/* MetaplotRegion functions */

MetaplotRegion::MetaplotRegion(int len)
{
	length = len;

	basePairs = new double*[len];
	
	for (int i = 0; i < len; i++)
	{
		basePairs[i] = new double[2]; // 0 = height, 1 = count
		basePairs[i][0] = 0;
		basePairs[i][1] = 0;
	}
}

void MetaplotRegion::addSignal(int offset, int len, double value, char strand)
{
	if (strand == '+')
	{
		for (int i = offset; i < len + offset; i++)
		{
		//	cerr << "DEBUG adding " << value << " to " << i << endl;
			basePairs[i][0] += value;
			basePairs[i][1]++;
		//	cerr << "basePairs " << i << " height is now " << basePairs[i][1] << " signal is " << basePairs[i][0] <<  endl;
		}		
	}
	else 
	{
		for (int i = length - offset - 1; i > length - (offset + len) - 1; i--)
		{
			basePairs[i][0] += value;
			basePairs[i][1]++;
		}
	}
}

/* End MetaplotRegion functions */

void calculateMetaplot(UserOpts &input, vector<Bed *> * &bedByChrs, vector<Wig *> &wigByChrs, vector<string> commonChrs, MetaplotRegion * &region)
{
//	cerr << "DEBUG in calculateMetaplot" << endl;
	for (vector<string>::iterator iter = commonChrs.begin(); iter != commonChrs.end(); iter++)
	{
		cout << endl << "On chromosome " << *iter << endl;
		// read in wig
		// then read in each bed

		cout << "Reading in wig file ... " << endl;
		vector<Wig *>::iterator wigIter;

		for (wigIter = wigByChrs.begin(); wigIter != wigByChrs.end(); wigIter++)
		{
			if ((*wigIter)->getChr() == (*iter))
				break;
		}
		(*wigIter)->readPeaks();

		cout << "\tdone." << endl;
		for (int bedFileNum = 0; bedFileNum < input.getBedNumber(); bedFileNum++)
		{
			cout << "\tOn bedfile " << bedFileNum << endl;

			cout << "\t\tReading in bedfile ... " << endl;
			vector<Bed *>::iterator bedIter;
			for (bedIter = bedByChrs[bedFileNum].begin(); bedIter !=bedByChrs[bedFileNum].end(); bedIter++)
			{
				if ((*bedIter)->getChr() == (*iter))
					break;
			}
			
			(*bedIter)->readPeaks();

			cout << "\t\t\tdone." << endl;

			Peak * currBedPeak = (*bedIter)->getCurrPeak();
			Peak * currWigPeak = (*wigIter)->getCurrPeak();
			
//			cerr << "DEBUG: segfault?" << endl;
			/* now begin actual processing */
			// assumes all peaks are in ASCENDING order (check the stacks!)
			// okay first hurdle: we need to go through the wig multiple times, but each bed only once. But they're all stacks.
			// unstack/restack wig
			 
			// call restack wig at the end of every bedfile processed
			// do NOT nextPeak() on wig file

			while(currBedPeak != NULL && currWigPeak != NULL)
			{
				// this loop runs through all bed peaks for this file/chr

				// find first wig peak that is not < bed peak
				while (currWigPeak != NULL && currBedPeak->start > currWigPeak->end)
				{
					(*wigIter)->unstack();
					currWigPeak = (*wigIter)->getCurrPeak();
				}
				// assumes wig peaks are smaller than bed peaks
				while(currWigPeak != NULL && currWigPeak->start < currBedPeak->end)
				{
					// match and add up
					int offset = 0;
					int len = 0;
					if (currWigPeak->start > currBedPeak->start)
					{
						offset = currWigPeak->start - currBedPeak->start;

						if (currWigPeak->end >= currBedPeak->end)
							len = currBedPeak->end - currWigPeak->start;
						else
							len = currWigPeak->end - currWigPeak->start;	
					}					
					else
					{
						offset = 0;
						
						if (currWigPeak->end > currBedPeak->end)
							len = currBedPeak->end - currBedPeak->start;
						else
							len = currWigPeak->end - currBedPeak->start;
					}
					region[bedFileNum].addSignal(offset, len, currWigPeak->value, currBedPeak->strand);
//					cerr << "DEBUG region[" << bedFileNum << "] has signal " << region[bedFileNum].basePairs[0][0] << " and count " << region[bedFileNum].basePairs[0][1] << endl;

					(*wigIter)->unstack();
					currWigPeak=(*wigIter)->getCurrPeak();
				}

				(*bedIter)->nextPeak();
				currBedPeak = (*bedIter)->getCurrPeak();
			}
			// now at end for this bedfile
			(*wigIter)->restack();
			(*bedIter)->clearPeaks();
		}
		(*wigIter)->clearPeaks();		
	}
//	cerr << "DEBUG end of calculateMetaplot"<<endl;
//	cerr << "DEBUG region[0][1] is " << region[0].basePairs[0][1] << endl;
}


void monteCarloMetaplot(string file, int bedNum)
{
 // avg horizontally
 // no need for int reps 	
	cerr << "DEBUG: in MonteCarloMetaplot" << endl;
	ifstream infile(file.c_str());
	ofstream outfile("metaplot_outfile.tmp");	

	outfile << "bp\tsimulation" << endl;

	string line;

	cerr << "DEBUG: before while" << endl;
	while (getline(infile, line))
	{
		// have to deal with NA's
		double num;
		stringstream linestream(line);

		cerr << "linestream is " << linestream.str() << endl;		

		if (!(linestream >> num))
			continue; // like next?
		
		// num == bp# right now
		double bp = num;
		double avg = 0;
		string NA;

// want to treat NAs like 0s
		for (int i = 0; i < bedNum; i++)
		{
			linestream >> NA;
			stringstream toNum(NA);

			if ((toNum >> num))
				avg += num;
		} 

		avg /= (double)bedNum;
		
		outfile << bp << "\t" << avg << endl;
	}
	cerr << "DEBUG: after while" << endl;
	
	outfile.close();
	infile.close();
	rename("metaplot_outfile.txt", "tmp");
	rename("metaplot_outfile.tmp", "metaplot_outfile.txt");
	rename("tmp", "metaplot_outfile.tmp");
	cerr << "DEBUG: end of function" << endl;
}


void debug(MetaplotRegion * &region, int bedNumber, string nameStr, string nameStrR)
{
	cerr << "In DEBUG function" << endl;
}

// returns outfile name
string printResults(MetaplotRegion * &region, int bedNumber, string nameStr, string nameStrR, int maxWindow)
{
//	debug();
	cout << endl << "Printing out results" << endl;

	string outfileName = "metaplot_outfile.txt";

	ofstream fstream(outfileName.c_str()); 

	string header = "bp\t" + nameStr + "\n";
	fstream << header;

	int h = 0 - (maxWindow / 2);

	for (int i = 0; i < maxWindow; i++)
	{
		int addition = h + i;
		stringstream hToString;
		hToString << addition;
		fstream << hToString.str() << "\t";

		for (int j = 0; j < bedNumber; j++)
		{
			if (region[j].basePairs[i][1] <= 0)
				fstream << "NA\t";
			else
			{
				double avg = (double)region[j].basePairs[i][0] / (double)region[j].basePairs[i][1];
				stringstream avgToString;
				avgToString << avg;
				fstream << avgToString.str() << "\t";
			}	
		}
		fstream << endl;
	}

	fstream.close();

	ofstream rstream("metaplot_outfile.R");
	rstream << "library(ggplot2)\nlibrary(reshape2)\n";
	rstream << "pdf(file=\"metaplot_outfile.pdf\", width=12, height=8)\n";
	rstream << "plot<-read.table(\"metaplot_outfile.txt\", header=T)\n";
	rstream << "plot.melt<-melt(plot[,c('bp', " << nameStrR << ")], id.vars=1)\n";
	rstream << "ggplot(plot.melt, aes(x=bp, y=value, colour=variable, group=variable)) +\n";
	rstream << "geom_line() +\n";
	rstream << "geom_smooth() +\n";
	rstream << "theme_bw() +\n";
	rstream << "ggtitle(\"Metaplot\") +\n";
	rstream << "theme(panel.grid.minor=element_blank()) +\n";
	rstream << "scale_colour_brewer(palette=\"Set1\", name=\"Bed\")\n";

	rstream.close(); 

	return outfileName;
}
