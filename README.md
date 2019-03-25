# Leri Analytics
[Leri Analytics](https://godzilla.uchicago.edu/pages/ngaam/leri/index.html) aims to provide bioinformatics services for and solutions to both academic labs and industries. Leri delivers better custom solutions in the fileds of 
- *in silico* protein folding & design (improved from [Sibe](https://godzilla.uchicago.edu/pages/ngaam/sibe/index.html))
- NGS data analysis
- population genomics
- statistical analysis on data
- deep learning & optimization
- as well as other omic data analysis.

![Alt text](leri-overview.png?raw=true "")

Leri Analytics delivers the most accurate solutions to your research and industrial applications. Understanding the specific purpose of your project, we will try our best to provide standard pipelines that particularly developed for your project, the best solutions, and the detailed analytical reports. Contact us now for a free consultation by email at yaan.jang@gmail.com.


## Build Leri from source
### 1. Dependencies
#### 1.1 Boost C++ libraries
```
sudo apt-get install libboost-all-dev
```
or
```
wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz 
tar -xvf boost_1_65_1.tar.gz 
cd boost_1_65_1 
./bootstrap.sh --libdir=/usr/local/lib --includedir=/usr/local/include 
./b2 
./b2 install 
```
#### 1.2 GFlags
```
wget https://github.com/schuhschuh/gflags/archive/master.zip 
unzip master.zip 
cd gflags-master 
mkdir build && cd build 
export CXXFLAGS="-fPIC" && cmake .. && make VERBOSE=1 
make && make install
```
#### 1.3 Glogs
```
wget https://github.com/google/glog/archive/v0.3.4.tar.gz 
tar zvxf glog-0.3.4.tar.gz 
cd glog-0.3.4/ 
./configure 
make && make install
```

### Build Leri
```
mkdir build
cd build
cmake ..
make
```
### APIs (part)
##### Base call
```cpp
int base_call() {
  printf("-- Leri Analytics: Base Calling\n");
  if (!FLAGS_fastx.size()) {
    printf("\033[1;31m");
    printf("ERROR: Where to locate your FAST5 data file?\n");
    printf("Usage:\n");
    printf("leri base_call -fastx <STRING> -type <INT> -output <STRING>");
    printf("\033[0m\n");
    exit(EXIT_FAILURE);
  }
  if (FLAGS_type < -1 && FLAGS_type > 3) {
    printf("\033[1;31m");
    printf("ERROR: What type of your data?\n");
    printf("Usage:\n");
    printf("-type=0 is Nanopore base calling\n");
    printf("-type=1 is Illumina base calling\n");
    printf("-type=2 is PacBio base calling\n");
    printf("leri base_call -fastx <STRING> -type <INT> -output <STRING>");
    printf("\033[0m\n");
    exit(EXIT_FAILURE);
  }
  printf("-- Invocation: leri base_call -fastx %s -type %i\n", FLAGS_fastx.c_str(), FLAGS_type);
  string atgc_name;
  string postfix = "leri.fastq";
  leri::define_output_file_name(atgc_name, FLAGS_fastx, FLAGS_output, postfix);
  // LOG(INFO)<<atgc_name;
  if (FLAGS_type == 0) {
    Nanopore<double> Fast5;
    Fast5.read_fast5(FLAGS_fastx, atgc_name);
  }
  if (FLAGS_type == 1) {  // TODO(Illumina):
    Fast5.read_fast5(FLAGS_fastx, atgc_name);
  }
  if (FLAGS_type == 2) {  // TODO(Pacbio):
    Fast5.read_fast5(FLAGS_fastx, atgc_name);
  }
  printf("-- Leri outputs: %s\n", atgc_name.c_str());
  return 0;
}
```
##### ATGC statistics
```cpp
/*===== ATGC =====*/
int atgc_stats() {
  if (!FLAGS_fastx.size()) {
    printf("\033[1;31m");
    printf("ERROR: Where to locate your FAST5 data file?\n");
    printf("Usage:\n");
    printf("leri atgc_stats -fastx <STRING> -type <INT> -output <STRING>");
    printf("\033[0m\n");
    exit(EXIT_FAILURE);
  }
  printf("-- Invocation: leri atgc_stats -fastx %s\n", FLAGS_fastx.c_str());
  ATGC<double> atgc(FLAGS_output, FLAGS_jobname);
  atgc.init();
  atgc.read_fastx(FLAGS_fastx);
  atgc.stats_report();
  return 0;
}
```
##### Read filtering
```cpp
int read_filter() {
  if (!FLAGS_fastx.size()) {
    printf("\033[1;31m");
    printf("ERROR: Where to locate your FASTQ data file?\n");
    printf("Usage:\n");
    printf("leri read_filter -fastx <STRING> -n <INT> (default 0) -output <STRING>\n");
    printf("Note: n is the allowed number of ambiguous characters\n");
    printf("\033[0m\n");
    exit(EXIT_FAILURE);
  }
  printf("-- Invocation: leri read_filter -fastx %s -n %i\n", FLAGS_fastx.c_str(), FLAGS_n);
  ATGC<double> atgc(FLAGS_output, FLAGS_jobname);
  atgc.init();

  string atgc_name;
  string postfix;
  postfix = "leri_filtered.fastq.gz";
  leri::define_output_file_name(atgc_name, FLAGS_fastx, FLAGS_output, postfix);
  atgc.read_filter(FLAGS_fastx, atgc_name, FLAGS_n);
  postfix = "leri_quality.fastq.gz";
  leri::define_output_file_name(atgc_name, FLAGS_fastx, FLAGS_output, postfix);
  atgc.read_quality(FLAGS_fastx, atgc_name, FLAGS_n);
  postfix = "leri_trimed.fastq.gz";
  leri::define_output_file_name(atgc_name, FLAGS_fastx, FLAGS_output, postfix);
  atgc.read_trim(FLAGS_fastx, atgc_name);
  return 0;
}
```
##### Kmer
```cpp
int kmer_generator() {
  if (!FLAGS_fastx.size()) {
    printf("\033[1;31m");
    printf("ERROR: Where to locate your FASTQ data file?\n");
    printf("Usage:\n");
    printf("leri kmer_generator -fastx <STRING> -min <INT> -max <INT> -output <STRING>\n");
    printf("\033[0m\n");
    exit(EXIT_FAILURE);
  }
  printf("-- Invocation: leri kmer_generator -fastx %s -min %i -max %i -output %s\n", FLAGS_fastx.c_str(), (int)FLAGS_min, (int)FLAGS_max, FLAGS_output.c_str());
  Kmer<double> kmer(FLAGS_fastx, FLAGS_output);
  kmer.generate_kmer_frequency((int)FLAGS_min, (int)FLAGS_max);
  return 0;
}
```
```cpp
int kmer_counter() {
  if (!FLAGS_fastx.size()) {
    printf("\033[1;31m");
    printf("ERROR: Where to locate your FASTQ data file?\n");
    printf("Usage:\n");
    printf("leri kmer_counter -type <INT> -fastx <STRING> -output <STRING>\n");
    printf("\033[0m\n");
    exit(EXIT_FAILURE);
  }
  printf("-- Invocation: leri kmer_counter -fastx %s\n", FLAGS_fastx.c_str());
  ATGC<double> fastq(FLAGS_output, FLAGS_jobname);
  fastq.init();
  fastq.kmer_counter(FLAGS_fastx);
  return 0;
}
```
##### Sequence assembly
```cpp
int atgc_assembly() {
  if (!FLAGS_fastx.size()) {
    info_error("ERROR: Where to locate your contig file?\nUsage: leri atgc_assembly -type <INT> -fastx <STRING> -output <STRING>\nwhere -type=0 is that Leri will use the De Bruijn Graph approach for assembly\n      -type=1 is that Leri will use the ACO algorithm for assembly\n");
  }
  printf("-- Invocation: leri atgc_assembly -type %d -fastx %s -kmer %i -output %s\n", FLAGS_type, FLAGS_fastx.c_str(), FLAGS_kmer,  FLAGS_output.c_str());
  printf("-- If you encounter an error of memery, try 'ulimit -u unlimited' in your terminal!\n");
  // Write reconstructed FASTQ into a file.
  string read;
  string atgc_name;
  string postfix = "leri_assembly.fastq";
  leri::define_output_file_name(atgc_name, FLAGS_fastx, FLAGS_output, postfix);
  leri::write_string(atgc_name, VERS_INFO("Read Assembly"), "w");
  if (FLAGS_type == 0) { // Assembly Based on de Bruijn Graph
    DeBrujinGraph<double> graph(FLAGS_fastx, FLAGS_kmer);
    graph.build_graph();
    // multimap<DeBrujinNode, DeBrujinNode>::iterator i;
    // for(auto i = gra.G.begin(); i != gra.G.end(); i++)
    //     cout << i->first.getKmer() << "-----" << i->second.getKmer() << endl;
    graph.assembled_sequence(read); 
  }
  if (FLAGS_type == 1) {  // Assembly Based on ACO algorithm
    for (int i = 0; i < FLAGS_ntrial; i++) {
      DNA<double> dna;
      dna.read_file(FLAGS_fastx);
      dna.init();
      dna.find_overlaps();
      vector<vector<int> > overlaps = dna.get_overlaps();
      vector<vector<int> > adj_mtx;
      int n_node = sqrt(overlaps.size()) + 1;
      //  num_ants, num_nodes, start_node, max_gen, ntrial
      vector<int> config = { FLAGS_npop, n_node, 0, FLAGS_max_iter };
      //  alpha, beta, rho, tau
      vector<double> param = { 0.8, 1.1, 0.15, 3.0, 4700.0 };

      leri::ACO<double> ants;
      ants.init(config, param, overlaps, adj_mtx);
      ants.actuation();
      vector<int> path = ants.get_best_path();
      dna.set_seq_order(path);
      string merged_seq = dna.join_fragments(adj_mtx);
      // printf("-- #%2i reconstruction: %s\n", i, merged_seq.c_str());
      merged_seq = std::to_string(i+1) + "  " + merged_seq;
      leri::write_string(atgc_name, merged_seq, "a+");
    }
  }

  leri::write_string(atgc_name, read, "a+");
  printf("-- Leri outputs: %s\n", atgc_name.c_str());
  return 0;
}
```
##### Statistical Tools
###### 1. Principle Component Analysis
```cpp
int stats_pca() {
  if (!FLAGS_mat.size()) {
    info_error("ERROR: Where to locate your data matrix file that requires PCA calculation?\nUsage: leri stat_pca -mat <STRING> -output <STRING>");
  }
  printf("-- Invocation: leri stats_pca -mat %s -output %s\n", FLAGS_mat.c_str(), FLAGS_output.c_str());
  bool all_normlzed = false;
  bool origin_normlzd = true;
  Stats_PCA<double> pca;
  pca.pca_init(FLAGS_mat);
  pca.pca_calc(all_normlzed, origin_normlzd);
  pca.pca_finish(FLAGS_mat, FLAGS_output);
  return 0;
}
```
###### 2. [Fast Independent Component Analysis](http://en.wikipedia.org/wiki/FastICA)
```cpp
int stats_ica() {
  if (!FLAGS_mat.size()) {
    info_error("ERROR: Where to locate your data matrix file that requires ICA calculation?\nUsage: leri stats_ica -mat <STRING> -output <STRING>");
  }
  printf("-- Invocation: leri stats_ica -mat %s -output %s\n", FLAGS_mat.c_str(), FLAGS_output.c_str());
  Stats_ICA<double> ica;
  ica.ica_calc(FLAGS_mat, FLAGS_output);
  ica.ica_finish(FLAGS_mat, FLAGS_output);
  return 0;
}
```
##### Sequentially evolving neural network (SENET)
```cpp
int senet() {
  printf( "-- Leri Analytics: Sequentially Evolving Neural Network on MNIST\n");
  printf( "-- Invocation: leri senet -threads %d -echo %d -mbs %d -solver %s -lrate %.3f -dirt %s\n", FLAGS_threads, FLAGS_echo, FLAGS_mbs, FLAGS_solver.c_str(), FLAGS_lrate, FLAGS_dirt.c_str());

  string str_tmp = FLAGS_dirt;
  int found = str_tmp.find_last_of("/\\");
  int len = str_tmp.length();
  if( found != len-1 ) { str_tmp += "/"; }
  FLAGS_dirt = str_tmp;

  leri::SENET<float> senet(FLAGS_dirt, FLAGS_output, FLAGS_jobname, FLAGS_solver, FLAGS_echo, FLAGS_mbs, FLAGS_lrate, FLAGS_threads);
  senet.senet_init();
  senet.senet_actuate();
  return 0;
}
```
##### OptiFel: fuzzy modeling method
```cpp
int optifel() {
  /* N. J. Cheung, X.-M. Ding, H.-B. Shen. IEEE Transactions on Fuzzy Systems, 22(4): 919-933, 2014. */
  FLAGS_type = 1;
  if (!FLAGS_mat.size()) {
    info_error("ERROR: Where to locate your data matrix file that requires for OptiFel modeling?\nUsage: leri optifel -mat <STRING> -max_iter <INT> -npop <INT> -nvar <INT> -nrule <INT> -output <STRING>");
  }
  if (FLAGS_phase != 0 && FLAGS_phase != 1) {
    info_error("ERROR: Leri doesn't know to train or predict, please tell it by -phase=0 (train) or -phase=1 (predict)\nUsage: leri optifel -mat <STRING> -max_iter <INT> -npop <INT> -nvar <INT> -nrule <INT> -output <STRING>");
  }
  if (!FLAGS_phase) {
    printf( "-- Invocation: leri optifel -mat %s -phase %d -max_iter %d -npop %d -nrule %d -output %s\n", FLAGS_mat.c_str(), FLAGS_phase, FLAGS_max_iter, FLAGS_npop, FLAGS_nrule, FLAGS_output.c_str());
  } else {
    printf( "-- Invocation: leri optifel -mat %s -cfg %s -phase %d -nrule %d -output %s\n", FLAGS_mat.c_str(), FLAGS_cfg.c_str(), FLAGS_phase, FLAGS_nrule, FLAGS_output.c_str());
  }

  leri::SolverParam<double> param;
  // param.m = m;
  param.set_param(FLAGS_nvar, FLAGS_npop);
  leri::CHPSO_Solver<double> OptiFel(FLAGS_type, param, FLAGS_output, FLAGS_jobname, FLAGS_ntrial, FLAGS_max_iter);
  if(FLAGS_phase == 0) { 
    OptiFel.optifel_load_data(FLAGS_mat, FLAGS_nrule);
    OptiFel.chpso_actuate();
  }
  if(FLAGS_phase == 1) {
    string str;
    leri::get_file_name(str, FLAGS_mat);
    OptiFel.optifel_load_data(FLAGS_mat, FLAGS_nrule, 1);
    OptiFel.optifel_load_data(FLAGS_cfg, FLAGS_nrule, 2);
    OptiFel.optifel_prediction(str, FLAGS_output);
  }
  OptiFel.finish();
  return 0;
}

```
##### Basic CNN classification
```cpp
int img_cnn() {
  printf( "-- Leri Analytics: Convolutional Neural Network on MNIST\n");
  if(!FLAGS_dirt.size()) {
    info_error("ERROR: Which directory to locate the image data! \nUsage: leri img_cnn -dirt <path-to-your-image-data>");
  }
  
  printf( "-- Invocation: leri img_cnn -echo %d -mbs %d -solver %s -lrate %.3f -dirt %s\n", FLAGS_echo, FLAGS_mbs, FLAGS_solver.c_str(), FLAGS_lrate, FLAGS_dirt.c_str());
  string str_tmp = FLAGS_dirt;
  int found = str_tmp.find_last_of("/\\");
  int len = str_tmp.length();
  if( found != len-1 ) { str_tmp += "/"; }
  FLAGS_dirt = str_tmp;

  leri::SENET<float> senet(FLAGS_dirt, FLAGS_output, FLAGS_jobname, FLAGS_solver, FLAGS_echo, FLAGS_mbs, FLAGS_lrate, FLAGS_threads);
  senet.senet_init();
  senet.senet_basic_cnn();
  return 0;
}  

```
##### Image format convertor
```cpp
int img2ubyte() {
  if (!FLAGS_dirt.size()) {
    printf("\033[1;31m");
    printf("ERROR: Which directory to locate the image data!\n");
    printf("Usage: leri img2ubyte -dirt <STRING> -w <INT> -h <INT> -gz <BOOL> -output <STRING>");
    printf("\033[0m\n");
    exit(EXIT_FAILURE);
  }
  printf("-- Invocation: leri img2ubyte -w %d -h %d -gz %d -dirt %s -output %s\n", FLAGS_w, FLAGS_h, FLAGS_gz, FLAGS_dirt.c_str(), FLAGS_output.c_str());
  string str_tmp = FLAGS_dirt;
  int found = str_tmp.find_last_of("/\\");
  int len = str_tmp.length();
  if (found != len-1) { str_tmp += "/"; }
  FLAGS_dirt = str_tmp;
  Image<double> im;
  im.convertor(FLAGS_dirt, FLAGS_output, FLAGS_w, FLAGS_h);
  return 0;
}

```
##### Protein Sequence Format Converter 
```cpp
int sequence_converter() {
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate the file (*.sto) generated by Jackhmmer?\n");
    printf( "Usage:\n");
    printf( "leri sequence_converter -msa=*' -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri sequence_converter -msa %s\n", FLAGS_msa.c_str()); 
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.jackhmmer2fasta(FLAGS_msa);
  sequence.finish();
  return 0;
}
```
##### Mix G and GPCR sequences
```cpp
int sequence_mix() {
  if(!FLAGS_fastx.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate G protein multiple sequence alignment (MSA) file?\n");
    printf( "Usage:\n");
    printf( "leri sequence_mix -fastx <G-protein-MSA> -msa <GPCR-MSA> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate GPCR multiple sequence alignment (MSA) file?\n");
    printf( "Usage:\n");
    printf( "leri sequence_mix -fastx <G-protein-MSA> -msa <GPCR-MSA> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf( "-- Invocation: leri sequence_mix -fastx %s -msa %s -output %s\n", FLAGS_fastx.c_str(), FLAGS_msa.c_str(), FLAGS_output.c_str());
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  
  sequence.initial();
  sequence.sequence_mix(FLAGS_fastx, FLAGS_msa);
  sequence.finish();
  
  return 0;
}
```
##### Protein Sequence Trimmer
```cpp
int sequence_trim() {
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate your multiple sequence alignment (MSA) file?\n");
    printf( "Usage:\n");
    printf( "leri sequence_trim -msa=* -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri sequence_trim -msa %s -gap_frac %.3f\n", FLAGS_msa.c_str(), FLAGS_gap_frac); 
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.read_msa(FLAGS_msa);
  sequence.sequence_trim(FLAGS_msa, FLAGS_gap_frac);
  sequence.finish();
  return 0;
}
```
##### Protein Sequence Statistics
```cpp
int sequence_stats() {
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate your multiple sequence alignment (MSA) file?\n");
    printf( "Usage:\n");
    printf( "leri sequence_stats -msa=*' -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri sequence_stats -msa %s\n", FLAGS_msa.c_str()); 
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.read_msa(FLAGS_msa);
  int nrow = sequence.get_num_of_seq();
  int nrow_scale = 2e3;
  if(nrow <= nrow_scale) {
    nrow_scale = nrow;
  } else {
    nrow_scale = 2e3;
    printf( "-- The size of the MSA is too large, but the similarities of the first %i sequences will be computed!\n", nrow_scale);
  }
  sequence.sequence_similarity(FLAGS_msa, nrow_scale);
  sequence.positional_conservation(FLAGS_msa);
  sequence.finish();
  return 0;
}
```
##### Residue evolutionary coupling
```cpp
int sequence_coupling() {
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to find your multiple sequence alignment (MSA) file?\n");
    printf( "Usage: leri sequence_coupling -msa <string> -output <string>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri sequence_coupling -msa %s\n", FLAGS_msa.c_str()); 

  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.read_msa(FLAGS_msa);
  #ifdef USE_GPU
  sequence.sequence_reweight_gpu();
  sequence.count_marginals_in_msa_gpu();
  sequence.pairwise_model_estimator_gpu();
  #else
  sequence.sequence_reweight();
  sequence.count_marginals_in_msa();
  sequence.pairwise_model_estimator();
  #endif

  sequence.write_potential();
  sequence.calc_dirct_couplings();
  sequence.calc_coupled_sector();
  sequence.finish();
  return 0;
}
```
##### Protein sequence energy
```cpp
int sequence_energy() {
  if(!FLAGS_msa.size() || !FLAGS_mat.size()) {
    printf( "\033[1;31m");
    if (!FLAGS_msa.size()) { 
      printf( "ERROR: where to locate your sequence(s) that require(s) energy calculation(s)?\n");
    }
    if (!FLAGS_mat.size()) { 
      printf( "ERROR: where to locate your sequence potential file?\n");
    }
    printf( "Usage:\n");
    printf( "leri sequence_energy -msa=* -mat=* -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri sequence_energy -msa %s -mat %s\n", FLAGS_msa.c_str(), FLAGS_mat.c_str()); 
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.read_msa(FLAGS_msa);
  sequence.calc_sequences_energy(FLAGS_msa, FLAGS_mat);
  sequence.finish();
  return 0;
}
```
##### Protein point mutation
```cpp
int point_mutation() {
  /* ddG = dG_mut-dG_wt */
  if(!FLAGS_fastx.size() || !FLAGS_mat.size()) {
    printf( "\033[1;31m");
    if (!FLAGS_fastx.size()) { 
      printf( "ERROR: where to locate your sequence that requires point mutaions?\n");
    }
    if (!FLAGS_mat.size()) { 
      printf( "ERROR: where to locate your sequence potential file?\n");
    }
    printf( "Usage:\n");
    printf( "leri point_mutation -fastx=* -mat=* -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri point_mutation -fastx %s -mat %s\n", FLAGS_fastx.c_str(), FLAGS_mat.c_str()); 
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.point_mutation(FLAGS_fastx, FLAGS_mat);
  sequence.finish();
  return 0;
}
```
##### Positionally conserved coupling
```cpp
int sequence_pcc() {
  printf("-- Leri Analytics: Positionally Conserved Coupling\n"); 
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to find your multiple sequence alignment (MSA) file?\n");
    printf( "Usage: leri sequence_pcc -msa <STRING> -threshold <FLOAT> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  if(FLAGS_threshold < 0) {
    printf( "\033[1;31m");
    printf( "ERROR: Please define a threshold to detect highly conserved positions. Default valude is 0.75.\n");
    printf( "Usage: leri sequence_pcc -msa <STRING> -threshold <FLOAT> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }

  printf("-- Invocation: leri sequence_pcc -msa %s -threshold %.3f\n", FLAGS_msa.c_str(), FLAGS_threshold); 
  string file_name;
  string postfix = "pcc";
  leri::define_output_file_name(file_name, FLAGS_msa, FLAGS_output, postfix);
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.calc_positionally_conserved_couplings(FLAGS_msa, file_name, FLAGS_threshold);
  sequence.finish();
  return 0;
}
```
##### GPCR-G barcode
```cpp
int gpcr_g_barcode() {
  printf("-- Leri Analytics: GPCR-G Protein Barcode\n"); 
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to find your multiple sequence alignment (MSA) file?\n");
    printf( "Usage: leri gpcr_g_barcode -msa <STRING>-output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri sequence_pcc -msa %s -gap_frac %.3f -cnsv_frac %.3f\n", FLAGS_msa.c_str(), FLAGS_gap_frac, FLAGS_cnsv_frac); 
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.read_msa(FLAGS_msa);
  sequence.sequence_trim(FLAGS_msa, FLAGS_gap_frac);
  sequence.sequence_positional_conservation_trim(FLAGS_msa, FLAGS_cnsv_frac);
  string aln_name, file_name;
  string postfix = "trimmed_pct.aln";
  leri::define_output_file_name(aln_name, FLAGS_msa, FLAGS_output, postfix);
  postfix = "positionally_conserved_coupling_sector_pct";;
  leri::define_output_file_name(file_name, FLAGS_msa, FLAGS_output, postfix);
  sequence.calc_positionally_conserved_couplings(aln_name, file_name);
  sequence.finish();
  return 0;
}
```
##### Protein sequence design
```cpp
int sequence_design() {
  if(!FLAGS_fastx.size() || !FLAGS_mat.size()) {
    printf( "\033[1;31m");
    if (!FLAGS_fastx.size()) {
      printf( "ERROR: where to locate your wild type sequence?\n");
    }
    if (!FLAGS_mat.size()) {
      printf( "ERROR: where to locate your sequence potential file?\n");
    }
    printf( "Usage:\n");
    printf( "leri sequence_design -fastx <STRING> -mat <STRING> -max_iter <INT>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri sequence_design -fastx %s -mat %s -max_iter %i, -temperature %.3f\n", 
             FLAGS_fastx.c_str(), FLAGS_mat.c_str(), FLAGS_max_iter, FLAGS_temperature); 

  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  int pos = FLAGS_fastx.find_last_of("/");
  printf( "-- WT sequence: %s, ", (FLAGS_fastx.substr(pos+1)).c_str());
  string str;// = FLAGS_fastx.substr(pos+1);
  leri::get_file_name(str, FLAGS_fastx);

  printf( "sequence trajactory: %s, ", (str+"_trajectory.txt").c_str());
  printf( "iterations: %i, ", FLAGS_max_iter);
  printf( "initial temperature: %.2f\n", FLAGS_temperature);
  sequence.set_max_iter(FLAGS_max_iter);
  sequence.set_temperature(FLAGS_temperature);
  sequence.sequence_design(FLAGS_fastx, FLAGS_mat, FLAGS_pdb);
  sequence.finish();
  return 0;
}
```
##### Estimate residue-contact
```cpp
int residue_contact() {
  if(!FLAGS_msa.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate your multiple sequence alignment (MSA) file?\n");
    printf( "Usage:\n");
    printf( "leri residue_contact -msa <STRING> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf("-- Invocation: leri residue_contact -msa %s\n", 
             FLAGS_msa.c_str()); 
  LSequence<double> sequence(FLAGS_output, FLAGS_jobname);
  sequence.initial();
  sequence.residue_contact(FLAGS_msa);
  sequence.finish();
  return 0;
}
```
##### PDB parser
```cpp
int pdb_parser() {
  if(!FLAGS_pdb.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: What is your PDB file?\n");
    printf( "Usage:\n");
    printf( "leri pdb_parser -pdb=* -type=* -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  if(!(FLAGS_type > -1 && FLAGS_type < 3)) {
    printf( "\033[1;31m");
    printf( "ERROR: What do you expect, FASTA or PBD or both of them?\n");
    printf( "Usage:\n");
    printf( "leri pdb_parser -pdb=* -type=* -output <STRING>\n");
    printf( "type=0 is to get fasta\ntype=1 is to reformat PDB\ntype=2 is to output both");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }

  int pos_1 = 0, pos_2 = 0;
  pos_1 = FLAGS_pdb.find_last_of("/");
  pos_2 = FLAGS_pdb.find_last_of(".");
  string tname = FLAGS_pdb.substr(pos_1+1);
  printf( "-- About the given PDB file: %s\n", tname.c_str());
  tname = FLAGS_pdb.substr(pos_1+1, pos_2-pos_1-1);

  Protein<double> protein;
  leri::ProteinData<double> model;
  model = protein.read_pdb(model, FLAGS_pdb);
  if (FLAGS_type == 0 || FLAGS_type == 2) {
    protein.write_fasta(FLAGS_output+tname+".fasta", model);
    printf( "-- Leri outputs FASTA file here: \n-- %s\n", (FLAGS_output+tname+".fasta").c_str());
  }
  if (FLAGS_type == 1 || FLAGS_type == 2) {
    protein.write_pdb(FLAGS_output+tname+"_leri.pdb", model);
    printf( "-- Leri outputs simplified PDB file here: \n-- %s\n", (FLAGS_output+tname+"_leri.pdb").c_str());
  }
  return 0;
}
```
##### Calculate residue contacts from a given PDB file
```cpp
int calc_contact() {
  if(!FLAGS_pdb.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: What is your PDB file?\n");
    printf( "Usage:\n");
    printf( "leri calc_contact -pdb=* -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  if (FLAGS_threshold == -1.0) {
    FLAGS_threshold = 7.5;
  }
  printf( "-- Invocation: leri calc_contact -threshold %.2f -atom_type %s -pdb %s\n", FLAGS_threshold, FLAGS_atom_type.c_str(), FLAGS_pdb.c_str());
  int pos_1 = 0, pos_2 = 0;
  pos_1 = FLAGS_pdb.find_last_of("/");
  pos_2 = FLAGS_pdb.find_last_of(".");
  string tname = FLAGS_pdb.substr(pos_1+1);
  printf( "-- About the given PDB file: %s\n", tname.c_str());
  tname = FLAGS_pdb.substr(pos_1+1, pos_2-pos_1-1);
  
  Protein<double> protein(FLAGS_output, FLAGS_jobname);
  leri::ProteinData<double> model;
  model = protein.read_pdb(model, FLAGS_pdb);
  protein.calc_residue_contact(tname, model, FLAGS_threshold, FLAGS_atom_type);
  return 0;
}

```
##### Extracting protein features from FASTA
```cpp
int build_feature_data() { // something wrong in reading mat3 files, should use read_func in utils.h
  if(!FLAGS_dirt.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: Where is your directory to load files?\n");
    printf( "Usage:\n");
    printf( "leri build_feature_data -dirt <STRING> -index <STRING> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  if(!FLAGS_index.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: What is your directory to load files?\n");
    printf( "Usage:\n");
    printf( "leri build_feature_data -dirt <STRING> -index <STRING> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }

  string tpath = FLAGS_dirt;
  int found = tpath.find_last_of("/\\");
  int len = tpath.length();
  if (found != len-1) { tpath += "/"; }
  FLAGS_dirt = tpath;

  string lib_path = string(LIB_DIRT) + "/";
  string cull_idx = FLAGS_dirt + FLAGS_index;//"cullpdb_pc50_res1.8_R0.25_d171102_chains10609_ok-all";
  string pdb_idx_name = FLAGS_output + "pdb.index";
  leri::collect_pdb(FLAGS_dirt, cull_idx, FLAGS_flag, pdb_idx_name);
  leri::seq2dat(lib_path, FLAGS_dirt, pdb_idx_name);
  //  FLAGS_w = 7;
  leri::build_data<double>(FLAGS_dirt, pdb_idx_name, FLAGS_w, FLAGS_output, "train");
  return 0;
}
```
##### Protein folding
```cpp
/*----- Protein Folding -----*/
int protein_folding() {
  if(!FLAGS_fastx.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate your FASTA sequence file to launch protein folding?\n");
    printf( "Usage:\n");
    printf( "leri protein_folding -fastx <STRING> -cfg <STRING> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  if(!FLAGS_cfg.size()) {
    printf( "\033[1;31m");
    printf( "ERROR: where to locate your paramter file to launch protein folding?\n");
    printf( "Usage:\n");
    printf( "leri protein_folding -fastx <STRING> -cfg <STRING> -output <STRING>");
    printf( "\033[0m\n");
    exit (EXIT_FAILURE);
  }
  printf( "-- Invocation: leri protein_folding -fastx %s -cfg %s -output %s\n", FLAGS_fastx.c_str(), FLAGS_cfg.c_str(), FLAGS_output.c_str());
  
  Simulation<double> simulation;

  simulation.select_driver_by_index(0); // There is only 1 driver now
  simulation.run_simulation(FLAGS_fastx, FLAGS_cfg.c_str(), FLAGS_output); // in simulation.cpp
  simulation.finish();
  return 0;
}
```
##### Optimization
###### 1. Convergent Heterogeneous Particle Swarm Optimizer
```cpp
int opt_chpso() {
  /* N. J. Cheung, X.-M. Ding, H.-B. Shen. IEEE Transactions on Fuzzy Systems, 22(4): 919-933, 2014. */
  printf( "-- Invocation: leri opt_chpso -max_iter %d -npop %d -nvar %d\n", FLAGS_max_iter, FLAGS_npop, FLAGS_nvar);
  FLAGS_type = 0;
  leri::SolverParam<double> param;
  // param.m = m;
  param.set_param(FLAGS_nvar, FLAGS_npop);
  leri::CHPSO_Solver<double> CHPSO(FLAGS_type, param, FLAGS_output, FLAGS_jobname, FLAGS_ntrial, FLAGS_max_iter);
  CHPSO.chpso_actuate();

  return 0;
}
```
###### 2. Differential Evolution (DE) method with self-adapting control parameters
```cpp
int opt_jde() {
  /* For more details, pleaser refer to 
     J. Brest, S. Greiner, B. Boskovic, M. Mernik, V. Zumer. IEEE Trans. Evol. Comput., 10(6): 646–657, 2006.
   */ 
  leri::SolverParam<double> param;
  // // param.m = m;
  leri::JDE_Solver<double> jDE(param);

  vector<vector<double> > ss;
  ss.resize(2, vector<double>(FLAGS_nvar, 0.0));

  int run_max = FLAGS_ntrial;  // Numbers of trials
  printf( "-- Invocation: leri opt_jde -max_iter %d -npop %d -nvar %d\n", FLAGS_max_iter, FLAGS_npop, FLAGS_nvar);
  printf( "     Number of varialbes: %i\n", FLAGS_nvar);
  double gopt;       // Best-so-far
  double avg = 0.0;
  int seed = rand() % 100000;               // In the range 0 to 99999;
  printf( "--  The jDE algorithm is to solve the given problem\n");
  jDE.jde_search_space(ss, FLAGS_nvar);
  for (int run = 0; run < run_max; run++) {
    seed += run;
    jDE.init(FLAGS_nvar, FLAGS_npop, seed);
    jDE.minimize(&gopt, ss, FLAGS_nvar, FLAGS_max_iter);
    avg += gopt;
    printf( "      #trial is %3i, best-so-far is %.8e\n", run+1, gopt);
    jDE.finish();
  }
  avg /= run_max;
  printf( "--  jDE: the average best-so-far of %i trials is %e\n", run_max, avg);
  return 0;
}
```
###### 3. Nested sampling
```cpp
int nested_sampler() {
  printf( "-- Invocation: leri nested_sampler -max_iter %d -npop %d -nvar %d\n", FLAGS_max_iter, FLAGS_npop, FLAGS_nvar);
  leri::ToyModel<double> model;
  leri::NestedSampler<double> sampler(FLAGS_npop, FLAGS_max_iter, model); // Particle number, steps, obj
  sampler.set_rng_seed(time(0));

  // // Do some iterations
  sampler.run(100*100);
  return 0;
}
```

## Examples
### 1. Protein sequence analysis
First, you need to install [HMMER software](http://hmmer.org) and download UniRef database [here](https://www.uniprot.org/downloads) to preprare protein sequences. 
Before you start to launch Leri, a multiple sequence alignment (MSA) is required. 
Firstly, install software *jackhmmer* on your own machine if there is no the tool. 
Here, *jackhmmer* is used to prepare the MSA by search against the Uniref\*\* database. 
A simple script that illustrates how to apply the *jackhmmer* to search your query 
protein sequence is show as follows,  

```
jackhmmer \
--notextwi \                                    # Unlimit ASCII text output line width
 -A <USER_DEFINED_FILE.sto> \                   # Save the multiple alignment of hits to file
--tblout <USER_DEFINED_FILE_tbl.out> \          # Save parseable table of per-sequence hits to file
--domtblout <USER_DEFINED_FILE_domtbl.out> \    # Save parseable table of per-domain hits to file
 -E 0.01 \                                      # Report sequences <= this E-value threshold in output
--popen 0.25 \                                  # Gap open probability
--pextend 0.4 \                                 # Gap extend probability
--mxfile <DIRECT_TO_BLOSUM62.mat> \             # Read substitution score matrix from file
 <FASTA> \                                      # Query protein sequence
 <DIRECT_TO_Uniref*> \                          # Where to locate the Uniref50, Uniref90, or Uniref100 database
>/dev/null
```

When you get the multiple alignment of hits from the jackhmmer, it is time to launch 
the Leri to convert the file to the stardard FASTA format. The command line that 
converts *.sto file to FASTA file is presented here,
```
leri sequence_converter -jobname <JOB_NAME> -msa <NAME_OF_STO>.sto
```
Trim aligned sequences according to the query sequence. 
```
leri sequence_trim -jobname <JOB_NAME> -msa <NAME_OF_MSA>_msa.a2m
```
Basic statistics on the aligned sequences,
```
leri sequence_stats -jobname <JOB_NAME> -msa <NAME_OF_MSA>_msa_trimmed.aln
```
### 2. Protein rational design
This project aims to provide a computational method for protein design. Accordingly, we can design a ``super protein" of a protein family stabilizing in a large range of temperature.
#### 2.1 Evolutionary coupling
```
leri sequence_coupling -jobname <JOB_NAME> -msa <NAME_OF_MSA>_msa_trimmed.aln 
```
#### 2.2 Single mutation
Leri provides a tool, shown as following command, for the computer-aided design of mutant proteins with controlled evolutionary information. It evaluates the changes at signle site in stability of a protein without referring to its tertiary structure. Fig. X illustrates energy changes (\Delta E=Emut-Ewt) between wild type and mutant sequences when a single mutation occurs to the wild type sequence.

#### 2.3 Rational design
```
leri sequence_design -jobname <JOB_NAME> -fastx <WILD_SEQ_NAME>.fasta -mat <POTENTIAL_NAME>.txt 
```
#### 2.4 Sequence energy
```
leri sequence_energy -jobname <JOB_NAME> -msa <SEQ_NAME>.aln -mat <POTENTIAL_NAME>.txt 
```

### 3. Sequencing data analysis
#### 3.1 Quality control
#### 3.2 Sequence assembly
### 4. Structural bioinformatics
#### 4.1 Get FASTA from PDB
```
leri pdb_parser -jobname <JOB_NAME> -pdb <PDB_NAME>.pdb -type <INT>
```
#### 4.2 Compute residue contacts
```
leri calc_contact -jobname <JOB_NAME> -pdb <PDB_NAME>.pdb -threshold 7.50 -atom_type CA 
```
#### 4.3 Protein folding
```
leri protein_folding -jobname <JOB_NAME> -fastx <FASTA_NAME>.fasta -cfg <CFG_NAME>.cfg 
```
#### 4.4 Features from FASTA
First to download ncbi-blast-2.7.1 and psipred-4.0.2, as well as database (nr database), then 
```
leri build_feature_data -jobname <JOB_NAME> -dirt <DIRECTORY> -flag 0 -index cullpdb_pc50_res1.8_R0.25_d171102_chains10609
```
where <index> is culled PDB ID list, one can download it from [PISCES server](http://dunbrack.fccc.edu/pisces), and <dirt> is to contain the index of PDB and results. 
#### 
### 5. Optimization
#### 5.1 The convergent heterogeneous particle swarm optimizer
The CHPSO algorithm is a swarm-based optimization approach. In CHPSO approach, four cooperative sub-swarms that can share information from each other but maintain differences contribute to keep a good balance between the exploration and exploitation in optimization and make the particles capable of converging to stable points. Users can define their objective function in include/leri/objfunc.hpp.
```
leri opt_chpso -jobname <JOB_NAME> -npop 8 -ntrial 3
```
#### 5.2 The Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm with limited-memory
The L-BFGS method is an optimization algorithm in the family of quasi-Newton methods using a limited amount of computer memory. Users can define their objective function in include/leri/objfunc.hpp.
#### 5.3 Nested sampling
```
leri nested_sampler -jobname <JOB_NAME> -npop 10 -max_iter 200 
```
#### 5.4 jDE
```
leri opt_jde -jobname <JOB_NAME> -npop 32
```
### 6. Machine learning
#### 6.1 Sequential evolving neural network (SENET)
SENET is a general evolving framework for deep neural networks on different specific data-sets.
#### 6.2 Phsior
Phsior is a real-value predictor developed based on convolutional neural network for predicting the torsion angles.
#### 6.3 OptiFel
OptiFel is a data-driven method of accurate and reliable fuzzy systems, where the structure and parameters are encoded in an optimization framework of heterogeneous heuristic search algorithm. This methodology provides different search behaviors for finding the optimal structure and parameters suitable for the subspaces of the fuzzy models. The theoretical proof has demonstrated that the cooperative search can maintain a balance between exploration and exploitation to ensure the algorithm converges to stable points. OptiFel can be applied to biological systems and protein structure prediction.

##### 6.3.1 OptiFel training phase
```
leri optifel -jobname <JOB_NAME> -mat leri-example/optifel/data-train.txt -nrule 3 -phase 0 -ntrial 3
```
##### 6.3.2 OptiFel test phase
```
leri optifel -jobname <JOB_NAME> -mat leri-example/optifel/data-test.txt -cfg <TRAINED_PARAMETERS>.par -nrule 3 -phase 1
```

#### 6.4 Basic CNN
To test simple CNN, one can use following command to start, and the report in HTML will be presented in *leri-output* if you don't define the output path. 
```
leri img_cnn -jobname CNN_Test -dirt leri-example/img/MNIST/ -threads <INT, option, e.g. 4> -echo <INT, option, e.g. 10> -mbs <INT, option, e.g. 32> -solver <STRING, option,e.g. adam> -lrate <FLOAT, option, e.g. 0.04>
```
You can also define your own output path using *-output <PATH_TO_YOUR_DIRECTORY>*.

History:

01-29-2019
Version 2.4.2 Leri Analytics: fixed bugs in de Brujin graph.

01-12-2019
Version 2.4.1 Leri Analytics: HTML and GNUPLOT based analytical reports are successfullly implemented.

12-18-2018
Version 2.4.0 Leri Analytics: protocols of sequential evolving neural network were initialized.

11-09-2018
Version 2.3.1 Leri Analytics: synteny is added and for testing.

11-02-2018
Version 2.3.0 Leri Analytics: GPU-ebabled calculations for protein sequences. 

10-27-2018
Version 2.2.9 Leri Analytics: assign GPUs to MPI process on PCC for sequence analysis.

10-25-2018
Version 2.2.9 Leri Analytics: assign GPUs to MPI process on PCC for sequence analysis.

09-29-2018
Version 2.2.8 Leri Analytics: fixed a bug in calculating positionally conserved couplings.

09-27-2018
Version 2.2.7 Leri Analytics: fixed a bug in reading MSA.

09-15-2018
Version 2.2.7 Leri Analytics: prototype of evolving neural network.

05-20-2018
Version 2.2.5 Leri Analytics: targeted to provide anlytical computation on big data and help researchers solve various data-ralated scientific problems.
 
03-29-2018 
Version 2.2.3 Leri Analytics: fixed bugs in Phsior and improve its performance.

03-26-2018 
Version 2.2.5 Leri Analytics: added a convertor from images to ubyte.

03-18-2018 
Version 2.2.4 Leri Analytics: added html reports for both deep learning and statistics on sequencing data.

02-27-2018 
Version 2.2.3 Leri Analytics: added operations including reads trim, filter, and quality score.

01-27-2018 
Version 2.2.2 Leri Analytics: statistical reports (in html) using plotly.

01-25-2018 
Version 2.2.1 Leri Analytics: add a module for DNA sequence assembly by ant colony optimization.

01-12-2018 
Version 2.2.0 Leri Analytics: work on sequencing data analysis.

01-09-2018 
Version 2.1.8. Leri Analytics: Fixed a bug in reading FASTQ.

01-06-2018 
Version 2.1.7. Leri Analytics: Fixed a bug in PCA.

10-07-2017 
Version 2.1.5 Leri Analytics: a website is developed and on test.

11-26-2017
Version 2.1.3 Leri Analytics: differential Evolution (DE) method with self-adapting control parameters has been included in optimizers.

09-25-2017 
Version 2.1.2 Leri Analytics: its architecture is completed and successful to work on OSX and UNIX-like system.

09-15-2017
Version 2.0.6 Leri Analytics: two main modules are added, (1) PCA and ICA for matrix (computed from MSA) analysis. (2) Basic statistical analysis for FASTQ and convert an SFF file from the 454 genome sequencer to FASTQ including the sequences and quality scores.

08-03-2017 
Version 2.0.0 Leri Analytics: older versoin is aborded. Now, it is projected to analytically learn from data and concern on Life Science. Support both Unix/Linux and OSX.
