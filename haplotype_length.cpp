//This program takes an aligned fasta file,a SNP position, and number of samples as input. Spits out minimum ancestral haplotype length for an allele. also provides minor allele count for that polymorphic site.


#include<iostream>
#include<fstream>
#include<vector>
#include<string>


using namespace std;

string cns_caller(vector<string> seq, int no_of_seq);
char max_elem(int *tot_count);
vector<string> read_fa(ifstream & fin);
vector<string> subsample(vector<string> seq, int sample_size, int pos, string cns);
int haplo_rcaller(vector<string> seq, string cns, int pos); // this can be converted into a void return; output will go into a file
int haplo_lcaller(vector<string> seq, string cns, int pos);


	int main(int argc, char * argv[])
	{
		if(argc ==1)
		{cerr<<"Usage: "<<argv[0]<<" seq_file no_of_samples position"<<endl;
		exit(EXIT_FAILURE);
		}

	ifstream fin;
	fin.open(argv[1]);
	vector<string> seq; 
	vector<string> s_seq; 
	string cns;  // store the consensus
	int h_pos; // store haplotype lengths

	int sample_size = atoi(argv[2]);
        cout<<"no. of samples "<<sample_size<<endl;
	int pos = atoi(argv[3]);
	cout<<"position given "<<pos<<endl;

	seq = read_fa(fin); //read all fasta sequences
	cout<<"Total no. of samples"<<seq.size()<<endl;
	cns = cns_caller(seq,sample_size); // call consensus by plurality for all sequences
	s_seq = subsample(seq,sample_size,pos,cns); // extract the sequences based on the allele at the polymorphic site
	cns = cns_caller(s_seq,s_seq.size()); // call consensus by plurality for the subsample
	
	cout<<"No. of alternate alleles "<<s_seq.size()<<endl;
	h_pos = haplo_rcaller(s_seq,cns,pos);
	cout<<"Right hand ancestral haplotype length "<<h_pos<<" bp"<<endl;
	h_pos = haplo_lcaller(s_seq,cns,pos);
	cout<<"Left hand ancestral haplotype length "<<h_pos<<" bp"<<endl;

	return 0;
	}

//****************functions***************
//****************************************
vector<string> read_fa(ifstream & fin)
{
string str;
vector<string> seq;
	while (getline(fin,str))
	{
		if(str[0] != '>')
		{seq.push_back(str);
		}
	}
return seq;
}


string cns_caller(vector<string> seq, int no_of_seq)
{
string str,str1;

int tot_count[] = {0,0,0,0}; //pos 1 for 'A',pos 2 for 'T',pos 3 for 'G', pos 4 for 'C'

str = seq[0]; // allocate the memory to str to save time.

	for(int i=0;i<str.size();i++) // for each position along the sequences
	{
		
		for(int j=0;j<no_of_seq;j++) // for each sequence
		{ str1 = seq[j];
			if(str1[i] == 'A')
			tot_count[0] = tot_count[0]+1;
			if(str1[i] == 'T')
			tot_count[1] = tot_count[1]+1;
			if(str1[i] == 'G')
			tot_count[2]=tot_count[2]+1;
			if(str1[i] == 'C')
			tot_count[3] = tot_count[3]+1;
		}

	str[i] = max_elem(tot_count);

	tot_count[0] =0; //reset the counts
	tot_count[1] =0;
	tot_count[2] =0;
	tot_count[3]=0;	
	}
return  str;
}


char max_elem(int *tot_count)
{
int n = tot_count[0];
char c = 'A';
	
	
			if(tot_count[1] > n)
			{ n = tot_count[1];
			c = 'T';
			}
			if(tot_count[2] > n)
			{ n = tot_count[2];
			 c = 'G';
			}
			if(tot_count[3] > n)
			{ n = tot_count[3];
			 c= 'C';
			}
	
return c;
}


vector<string> subsample(vector<string> seq,int sample_size, int pos, string cns)  // returns a subsample of a string vector which does not have cns[i] at position i 
{
vector<string> subsequence;
string str;
	for(int i =0;i<sample_size;i++)
	{	
	str = seq[i];
		if((str[pos-1] != cns[pos-1]) && (str[pos-1] != 'N') && (( str[pos-1] =='A')||( str[pos-1]=='T')||(str[pos-1]=='G')||(str[pos-1]=='C'))) // pos-1 is used because vector index starts at 0 and not 1; don't use the sequences with ambiguous codes
		{subsequence.push_back(str);
		}
	}
return subsequence;
}


int haplo_rcaller(vector<string> seq, string cns, int pos)
{

int h_pos;
string str; // to hold a sequence from the vector
	for(int i = pos-1;i<cns.size();i++) // start from polymorphism and slide to the right
		{
		for(int j=0;j<seq.size();j++)
			{
			 str = seq[j];
				if((str[i] != cns[i])&&(str[i]!='N')) // if a mismatch with the consensus (i.e. polymorphism) of the subsample is found
				{h_pos = (i-pos+1);
				i = cns.size();
		
				}
			}
		}
return h_pos;
}


int haplo_lcaller(vector<string> seq, string cns, int pos)
{
int h_pos;
string str;
	for(int i=pos-1;i>=0;--i) // start from the polymorphism and slide to the left
	{
		for(int j=0;j<seq.size();j++)
		{
		str = seq[j];
			if(str[i] != cns[i] && (str[i] != 'N'))
			{
			h_pos = (pos-1 -i);
			i = -1;
			}
		}
	}
return h_pos;
}
