#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
/*
 * Read in list of fasta files [ before multi-sequence alignment ]
 * Find the Ns in the fasta files
 * Remove the Ns and concatenate these sequences
 * Output them into a similar fasta file 
 */

using namespace seqan;

template <typename TString>
String<Dna> remove_ns(TString const & seq)
{
 String<Dna> no_n_seq;
 resize(no_n_seq, length(seq));
 unsigned no_n_seq_pos = 0;
 for (unsigned i = 0; i < length(seq); ++i)
 {
  if (seq[i] != 'N' && seq[i] != 'n')
  {
    no_n_seq[no_n_seq_pos] = seq[i];
    ++no_n_seq_pos;
  }
 } 
 return no_n_seq; 
}



int main(int argc, char const **argv) {  
  ArgumentParser parser("remove_ns");//
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "TEXT"));//
  ArgumentParser::ParseResult res = parse(parser, argc, argv);//
  if (res != ArgumentParser::PARSE_OK)
    return res == ArgumentParser::PARSE_ERROR;
  
  String<char> filepath;
  getArgumentValue(filepath, parser, 0);//
  String<char> qName;
  String<Dna5> seq;
  String<Dna> no_n_seq;
  
  String<char> input_filepath = filepath;
  append(input_filepath, ".fa");
  String<char> output_filepath = filepath;
  append(output_filepath, "_rm_ns.fa");
  
  SequenceStream seqStream(toCString(input_filepath));//
  if (!isGood(seqStream))
  {
    std::cerr << "ERROR: Failed to open input file." << std::endl;
    return 1;
  }  
  SequenceStream outputSeqStream(toCString(output_filepath), SequenceStream::WRITE);//
  if (!isGood(outputSeqStream))
  {
    std::cerr << " ERROR: Failed to open output file." << std::endl;
    return 1;
  }
  
  while (!atEnd(seqStream))
  {    
   clear(no_n_seq); 
   clear(seq);
   if (readRecord(qName, seq, seqStream) != 0)//
   {
     std::cerr << "ERROR: Failed to read record: " << qName << std::endl;
     return 1;     
   }   
   no_n_seq = remove_ns(seq);   //
   if (writeRecord(outputSeqStream, qName, no_n_seq) != 0)//
   {
     std::cerr << " ERROR: Failed to write record: " << qName << std::endl;
     return 1;
   }      
  }  
  return 0;  
}
