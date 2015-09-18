#ifndef METAPLOT_H
#define METAPLOT_H

#include<iostream>
#include<string>
#include<stack>
#include<vector>
#include<iterator>
#include<exception>
#include<sstream>
#include<fstream>
#include<cstdio>
#include<dirent.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>


typedef struct Peak
{
	int start;
	int end;
	double value;
	char strand;
	std::string chr;
} Peak;

class Chromosome
{
	public:
	Chromosome() {};
	Chromosome(std::string chr) 
	{
		chrName = chr;
	}	

	virtual void addPeak(void);
	virtual void printPeaks(void); // delete peak stack!
	virtual void readPeaks(void);
	
	std::string getChr(void);
	void nextPeak(void);
	Peak * getCurrPeak(void);	
	int getPeakStart(void);
	int getPeakEnd(void);
	void setFilename(std::string filename);
	void clearPeaks(void);	

	protected:
	virtual void generateFilename(void);

	std::string chrName;
	std::string fileName;
	std::stack<Peak> peaks;
};

class Wig: public Chromosome
{
	public:
	Wig() {};
	Wig(std::string chr) : Chromosome(chr) { };
	virtual void printPeaks(void); // delete peak stack!
	virtual void readPeaks(void);
	virtual void addPeak(int start, int end, double value);
	virtual void generateFilename(void);

	void unstack(void);
	void restack(void);

	private:
	
	std::stack<Peak> unstacked;
};

class Bed : public Chromosome
{
	public:
	Bed() {};
	Bed(std::string chr, int len) : Chromosome(chr)
	{
		peakLen = len;
	}

	virtual void printPeaks(void); // delete peak stack!
	virtual void readPeaks(void);
	virtual void addPeak(int start, int end, char strand);
	virtual void generateFilename(int filenum);
	
	private: 
	int peakLen;
};

class MetaplotRegion
{
	public:
	MetaplotRegion() { };
	MetaplotRegion(int len);

	void resetIndicies(void);
	void addSignal(int offset, int len, double value, char strand);

	double ** basePairs;
	private:
	int length;
};

typedef struct Freq
{
	std::string chr;
	int freq;
} Freq;

class UserOpts
{
	public:
	UserOpts(int argc, char *argv[]);

	void printUsage(void);
	int getBedNumber(void);
	int getMaxWindow(void);
	std::string getNextBedFile(void);
	void resetIter(void);
	std::vector<Bed *> * splitBedFiles(void);
	std::vector<Wig *> splitWigFile(void);
	void debugVectorSize(void);

	std::string getNameString(void);
	std::string getNameStringR(void);
	std::vector<std::string> commonChroms(void);
	
	private:
	bool validateBeds(void);
	bool validateWig(void);
	bool validateNames(void);

	std::string handleOpts(int argc, char * argv[]);
	void preprocessBed(void);
	void preprocessWig(void);

	void calcSlidingWindow(std::vector<Peak> * wigBlocks);

	int numBeds;
	int maxWindow;

	int step;	
	int window;
	int mode;

	bool preprocessBedOpt;
	bool preprocessWigOpt;

	std::string splitWigDir;
	std::string splitBedDir;
	bool wigsSplit;	
	bool bedsSplit;

	std::vector<Bed *> * readBedSplit(std::vector<Bed *> * &arOfBedfiles);
	std::vector<Wig *> readWigSplit(std::vector<Wig *> &wigSplit);

	std::stack<std::string> getFilesInDir(std::string dir);
	
	std::vector<Freq> chroms;	

	std::string wigFile;
	std::vector<std::string> bedFiles;
	std::stack<std::string> names;
	std::stack<std::string> namesR;
	// add more like debug mode later
};

void calculateMetaplot(UserOpts input, std::vector<Bed *> * bedByChrs, std::vector<Wig *> wigByChrs, std::vector<std::string> commonChrs, MetaplotRegion * region);

#endif
