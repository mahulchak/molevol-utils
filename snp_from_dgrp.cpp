#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

using namespace std;

int report_count(string str,string substr);
int cord_transform(string str);

int main()

{
string str;
ifstream fin;
fin.open("dgrp2.vcf");
int seq_ct =0 ;
bool start=0;
size_t found,found1 =0;

cout<<"Your coordinates are (9385388 is pos 1)\n";

	while(getline(fin,str))
	{
	found = str.find("2L_9385390_SNP");
		 
		if(found == 11)
		{start =1;}
		
		if(start == 1)
		{
			if(report_count(str,"ALTCOUNT=") <40 && report_count(str,"ALTCOUNT=")>5 && str.find("DEL") == std::string::npos)
			{
			seq_ct++;
			//cout<<str<<endl;
			cout<<cord_transform(str)<<endl;
			}
		}
	found1 = str.find("2L_9393544_SNP") ;
                 if(found1 == 11)
                 {break;}
	
	}
cout<<" Total number of sequences "<<seq_ct<<endl;

fin.close();
return 0;
}



int report_count(string str, string substr)

{
size_t found,found1;
int n;

found = str.find(substr);
found1 = str.find("\t",found);
n = found+substr.size();

string num = str.substr(n,(found1-n)); 
n = atoi(num.c_str());

return n;
}

int cord_transform(string str)
{
int n = atoi(((str.substr(3,10)).c_str()));
n = (n-9385388); 
return n;
}
	
