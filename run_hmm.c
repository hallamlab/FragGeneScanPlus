#include "hmm.h"
#include "util_lib.h"
#include "fasta.h"
#include <stdio.h>

unsigned int ids; // thread id
unsigned int threadnum = 1;
unsigned int format;
unsigned int wholegenome;
unsigned int output_dna;
unsigned int output_meta;
unsigned int verbose;
unsigned int MAX_BYTES_PER_BUFFER;
unsigned int MAX_SEQS_PER_BUFFER;
unsigned int num_reads_flag = 0;

pthread_t writer_thread;
long long round_counter;
off_t stopped_at_fpos; // tracks how far we've read in the input 

long writer_counter = 0;
long read_counter = 0;
long read_counter1 = 0;
long work_counter = 0;
long total_reads = -1;
long viterbi_counter = 0;
long num_writes = 0, num_reads=0;

SEM_T work_sema;
SEM_T stop_sema;
SEM_T counter_sema;

FILE* outfile_fp;
FILE* dna_outfile_fp;

char mystring[STRINGLEN];
char complete_sequence[STRINGLEN];

FASTAFILE* fp;

void setupProgram(int argc, char** argv){

  int max_mem = 0;
  fp = 0;
  // strncpy(train_dir, argv[0], strlen(argv[0])-4);
  // Set default training directory
  strcpy(train_dir, "train");
  strcat(train_dir, "/");
  strcpy(mstate_file, train_dir);
  strcat(mstate_file, "gene");
  strcpy(rstate_file, train_dir);
  strcat(rstate_file, "rgene");
  strcpy(nstate_file, train_dir);
  strcat(nstate_file, "noncoding");
  strcpy(sstate_file, train_dir);
  strcat(sstate_file, "start");
  strcpy(pstate_file, train_dir);
  strcat(pstate_file, "stop");
  strcpy(s1state_file, train_dir);
  strcat(s1state_file, "stop1");
  strcpy(p1state_file, train_dir);
  strcat(p1state_file, "start1");
  strcpy(dstate_file, train_dir);
  strcat(dstate_file, "pwm");
  /* read command line argument */
  //!! This argument reading should all be encapsulated in a single function, this will make reading the code much easier, right now we have to always move around it.
  if (argc <= 8){    
    fprintf(stderr, "ERROR: You missed some parameters for input\n");
    print_usage();
    exit(EXIT_FAILURE);
  }

  int c;    

  while ((c=getopt(argc, argv, "fs:m:o:w:r:t:p:dev")) != -1){
    switch (c){
      case 's':
        strcpy(seq_file, optarg);
        if (access(seq_file, F_OK)==-1){
          fprintf(stderr, "ERROR: Sequence file [%s] does not exist\n", seq_file);
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;  
      case 'm':
        max_mem = atoi(optarg);
        break;
      case 'w':
        wholegenome = atoi(optarg);
        if (wholegenome != 0 && wholegenome != 1){
          fprintf(stderr, "ERROR: An incorrect value for the option -w was entered\n");
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;
      case 'p':
        threadnum = atoi(optarg);
        if (threadnum < 1){
          fprintf(stderr, "ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        strcpy(out_file, optarg);
        strcpy(aa_file, out_file);
        strcat(aa_file, ".faa");
        strcpy(dna_file, out_file);
        strcat(dna_file, ".ffn");
        break;
      case 'r':
        // Location of train directory
        strcpy(train_dir, optarg);
        strcat(train_dir, "/");
        strcpy(mstate_file, train_dir);
        strcat(mstate_file, "gene");
        strcpy(rstate_file, train_dir);
        strcat(rstate_file, "rgene");
        strcpy(nstate_file, train_dir);
        strcat(nstate_file, "noncoding");
        strcpy(sstate_file, train_dir);
        strcat(sstate_file, "start");
        strcpy(pstate_file, train_dir);
        strcat(pstate_file, "stop");
        strcpy(s1state_file, train_dir);
        strcat(s1state_file, "stop1");
        strcpy(p1state_file, train_dir);
        strcat(p1state_file, "start1");
        strcpy(dstate_file, train_dir);
        strcat(dstate_file, "pwm");
        break;
      case 't':
        strcpy(train_file, optarg);
        strcpy(hmm_file, train_dir);
        strcat(hmm_file, train_file);
        if (access(hmm_file, F_OK)==-1){
          fprintf(stderr, "ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;
      case 'f':
        format = 1;
        break;
      case 'd':
        output_dna = 1;
        break;
      case 'e':
        output_meta = 1;
        break;
      case 'v':
        verbose = 1;
        break;
    }
  }

  optind = 1;
  while ((c=getopt(argc, argv, "fs:m:o:w:r:t:p:dev")) != -1){
    switch (c){
      case 'r':
        // Location of train directory
        strcpy(train_dir, optarg);
        strcat(train_dir, "/");
        strcpy(mstate_file, train_dir);
        strcat(mstate_file, "gene");
        strcpy(rstate_file, train_dir);
        strcat(rstate_file, "rgene");
        strcpy(nstate_file, train_dir);
        strcat(nstate_file, "noncoding");
        strcpy(sstate_file, train_dir);
        strcat(sstate_file, "start");
        strcpy(pstate_file, train_dir);
        strcat(pstate_file, "stop");
        strcpy(s1state_file, train_dir);
        strcat(s1state_file, "stop1");
        strcpy(p1state_file, train_dir);
        strcat(p1state_file, "start1");
        strcpy(dstate_file, train_dir);
        strcat(dstate_file, "pwm");
        break;
    }
  }


  /* check whether the specified files exist */
  if (access(mstate_file, F_OK)==-1){
    fprintf(stderr, "Forward prob. file [%s] does not exist\n", mstate_file);
    exit(1);
  }
  if (access(rstate_file, F_OK)==-1){
    fprintf(stderr, "Backward prob. file [%s] does not exist\n", rstate_file);
    exit(1);
  }
  if (access(nstate_file, F_OK)==-1){
    fprintf(stderr, "noncoding prob. file [%s] does not exist\n", nstate_file);
    exit(1);
  }
  if (access(sstate_file, F_OK)==-1){
    fprintf(stderr, "start prob. file [%s] does not exist\n", sstate_file);
    exit(1);
  }
  if (access(pstate_file, F_OK)==-1){
    fprintf(stderr, "stop prob. file [%s] does not exist\n", pstate_file);
    exit(1);
  }
  if (access(s1state_file, F_OK)==-1){
    fprintf(stderr, "start1 prob. file [%s] does not exist\n", s1state_file);
    exit(1);
  }
  if (access(p1state_file, F_OK)==-1){
    fprintf(stderr, "stop1 prob. file [%s] does not exist\n", p1state_file);
    exit(1);
  }
  if (access(dstate_file, F_OK)==-1){
    fprintf(stderr, "pwm dist. file [%s] does not exist\n", dstate_file);
    exit(1);
  }
  if (access(hmm_file, F_OK)==-1){
    fprintf(stderr, "hmm file [%s] does not exist\n", hmm_file);
    exit(1);
  }

  /* check for mem limit, allocate buffer */
  if (max_mem <= 0) {
    printf("Max memory limit specified invalid, defaulting to 1024MB\n");
    max_mem = 1024;
  } 

  // 5 stands for the number of buffers we are currently using per thread
  MAX_BYTES_PER_BUFFER = max_mem*100000/(5*2*threadnum);
  MAX_SEQS_PER_BUFFER = MAX_BYTES_PER_BUFFER/STRINGLEN;

  // remove them, if they already exist
  remove(aa_file);
  if (output_meta) remove(out_file);
  if (output_dna) remove(dna_file);
}

void initializeSemaphores(){

#ifdef __APPLE__
  sem_unlink("/work_sema");
  if ( ( work_sema = sem_open("/work_sema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

  sem_unlink("/sema_Q");
  if (( sema_Q = sem_open("/sema_Q", O_CREAT, 0644, 1))== SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

  sem_unlink("/sema_R");
  if (( sema_R = sem_open("/sema_R", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

  sem_unlink("/sema_r");
  if (( sema_r = sem_open("/sema_r", O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

  sem_unlink("/sema_w");
  if (( sema_w = sem_open("/sema_w", O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

  sem_unlink("/stop_sema");
  if (( stop_sema = sem_open("/stop_sema", O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

  sem_unlink("/COUNTER_SEMA");
  if (( counter_sema = sem_open("/COUNTER_SEMA", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

#elif __linux
  sem_init(&work_sema, 0, 1);
  sem_init(&sema_Q, 0, 1);   
  sem_init(&sema_R, 0, 1);
  sem_init(&sema_r, 0, 0);
  sem_init(&sema_w, 0, 0);
  sem_init(&stop_sema, 0, 0);
  sem_init(&counter_sema, 0, 1);
#endif
}

void destroySemaphores(){

#ifdef __APPLE__
  sem_unlink("/work_sema");
  sem_unlink("/sema_Q");
  sem_unlink("/sema_R");
  sem_unlink("/sema_r");
  sem_unlink("/sema_w");
  sem_unlink("/stop_sema");
  sem_unlink("/COUNTER_SEMA");

  char name[40];
  int j;
  for(j=0; j<threadnum; j++) {
    sprintf(name, "/sema_r%d", j);
    sem_unlink(name);

    sprintf(name, "/sema_w%d", j);
    sem_unlink(name);
  }
#endif

}

void initializeThreads(){

  pthread_t *thread = malloc(sizeof(pthread_t*) * threadnum);
  memset(thread, 0, sizeof(pthread_t*) * threadnum);
  thread_datas = malloc(sizeof(thread_data) * threadnum);

  // allocate memory for each thread only once!
  if (verbose)
    printf("DEBUG: Allocating memory for all threads...\n");

  int k;
  for(k=0; k< threadnum; k++) {
    init_thread_data(thread_datas+k);
  }

  if (verbose) {
    printf("DEBUG: Allocated memory for all threads!\n");
    printf("DEBUG: Starting writer thread...\n");
  }

  pthread_create(&writer_thread, 0, writer_func, 0);

  fp = OpenFASTA(seq_file);

  if (!fp) {
    printf("ERROR! Could not open seq_file %s for reading...!\n", seq_file);
    exit(0);
  }

  if (verbose)
    printf("INFO : Giving workers initial inputs...\n");

  int i,j;
  void *status;
  for(j=0; j<threadnum; j++)
    pthread_create(&thread[j], 0, thread_func, (void *)(thread_datas+j));

  for(j=0; j<threadnum; j++) {
    for(i=0;i<2;i++) {
      if( (stopped_at_fpos = read_seq_into_buffer(fp, thread_datas + j, i))!=0){
        sem_post(thread_datas[j].sema_r);                    
      }
    }
  }

  if (verbose)
    printf("INFO : Initializing worker threads...\n");

}

void conductWork(){
  // master loop - while we haven't exhausted reading the file yet
  while (stopped_at_fpos!=0) {
    sem_wait(sema_r);

    sem_wait(sema_Q);        
    QUEUE* temp;
    cutnpaste_q(&temp, EMPTY_Q);
    sem_post(sema_Q);

    while(temp) {
      sem_wait(sema_R);
      stopped_at_fpos = read_seq_into_buffer(fp,  temp->td, temp->buffer);
      sem_post(sema_R);

      sem_post(temp->td->sema_r);
      temp = temp->next;
    }
  }
  CloseFASTA(fp);

  if (verbose)
    printf("INFO : Finished handing out all the work...\n");

  num_reads_flag =1;

  sem_wait(stop_sema);
}

int main (int argc, char **argv) {

  // read in options and setup the necessary variables
  setupProgram(argc, argv);

  // setup all of the semaphores
  initializeSemaphores();

  if (verbose)
    printf("DEBUG: Max number of sequences per thread : %d, max bytes per thread : %d\n", MAX_SEQS_PER_BUFFER, MAX_BYTES_PER_BUFFER*5);

  /* read all initial model */
  get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

  // prepare all of the worker threads as well as the writer thread
  initializeThreads();

  // master loop - while we haven't exhausted reading the file yet
  conductWork();

  // destroy the semaphores if we have a mac machine
  destroySemaphores();

  printf("Run finished with %d threads.\n", threadnum);
}

int read_seq_into_buffer(FASTAFILE* ffp, thread_data* thread_data, unsigned int buf) {

  char *seq;
  char *name;
  int   L;
  int status ;
  int count=0;

  while ( (count < MAX_SEQS_PER_BUFFER) && (status = ReadFASTA(ffp, &seq, &name, &L)) ==1 ) {
    strcpy(thread_data->input_head_buffer[buf][count], name);
    strcpy(thread_data->input_buffer[buf][count], seq);
    read_counter++;
    count++;
  }

  thread_data->input_num_sequences[buf] = count;
  read_counter1 += thread_data->input_num_sequences[buf];

  return count;

}

// read as much of the file as you can into the buffer, return the file position
// of the last carat character encountered before the memory limit was hit
off_t read_file_into_buffer(FILE* fp, int fpos, thread_data* thread_data, unsigned int buf) {

  int i = 0; // temp length of each sequence
  long count = 0; // index for the number of sequences
  off_t last_carat_position = 0;
  long long total_count = 0;
  int seq_length; 

  memset(complete_sequence, 0, STRINGLEN);
  memset(mystring, 0, STRINGLEN);

  // loads the input
  while ( fgets(mystring, STRINGLEN, fp) && count < MAX_SEQS_PER_BUFFER-1 ) { 
    if (mystring[strlen(mystring)-1] == '\n') mystring[strlen(mystring)-1] = 0;
    if (mystring[0] == '>'){
      if (i>0){
        total_count += i;
        memset(thread_data->input_buffer[buf][count], 0, STRINGLEN);
        strcpy(thread_data->input_buffer[buf][count], complete_sequence);
        printf("%s\n",thread_data->input_buffer[buf][count]);
        memset(complete_sequence, 0, STRINGLEN);
      }   
      count++;
      read_counter++;
      strcpy(thread_data->input_head_buffer[buf][count], mystring);
      last_carat_position = ftello(fp);
      i = 0;
    } else {
      strcat(complete_sequence, mystring);
      seq_length = strlen(mystring);
      while(mystring[seq_length-1] == 10 || mystring[seq_length-1]==13 || mystring[seq_length-1]==0){
        seq_length --; 
      }   
      i += seq_length;
    }   
  }   

  strcpy(thread_data->input_buffer[buf][count], complete_sequence);
  thread_data->input_num_sequences[buf] = count;
  read_counter1 += thread_data->input_num_sequences[buf];

  if (ferror(fp)) {
    printf("Error with file encountered.\n");
  } else if(feof(fp)) {
    last_carat_position = -1; // terminating condition
  }   
  num_reads++;

  return last_carat_position;
}


void init_thread_data(thread_data* td){
  // Initialize thread data structure

  td->hmm = malloc(sizeof(HMM));
  memset(td->hmm, 0, sizeof(HMM));
  memcpy(td->hmm, &hmm, sizeof(HMM));

  td->wholegenome = wholegenome;

  td->id = ids;
  ids++;

#ifdef __APPLE__
  char name[40];
  sprintf(name, "/sema_r%d", td->id);
  sem_unlink(name);

  if((td->sema_r = sem_open(name, O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

  sprintf(name, "/sema_w%d", td->id);
  sem_unlink(name);
  if (( td->sema_w = sem_open(name, O_CREAT, 0644, 2)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }

#elif __linux
  sem_init(&td->sema_r, 0, 0);
  sem_init(&td->sema_w, 0, 2);
#endif

  // TODO : refactor to as many single large malloc calls as possible
  td->output_num_sequences = malloc(sizeof(int) * 2);
  memset(td->output_num_sequences, 0, sizeof(int)*2);
  td->input_num_sequences = malloc(sizeof(int) * 2);
  memset(td->input_num_sequences, 0, sizeof(int)*2);

  td->input_buffer = (char ***)malloc(sizeof(char**) * 2);
  td->input_head_buffer = (char ***)malloc(sizeof(char**) * 2);    
  td->output_buffer = (char ***)malloc(sizeof(char**) * 2);    
  td->aa_buffer = (char ***)malloc(sizeof(char**) * 2);    
  td->dna_buffer = (char ***)malloc(sizeof(char**) * 2);    
  td->acceptable_buffer = (unsigned int **) malloc(sizeof(unsigned int*) * 2);    

  td->dna	= malloc(sizeof(char) * 1500);
  memset(td->dna, 0, sizeof(char) * 1500);
  td->dna1 = malloc(sizeof(char) * 1500);
  memset(td->dna1, 0, sizeof(char)*1500);
  td->dna_f = malloc(sizeof(char) * 1500);
  memset(td->dna_f, 0, sizeof(char)*1500);
  td->dna_f1 = malloc(sizeof(char) * 1500);
  memset(td->dna_f1, 0, sizeof(char)*1500);
  td->protein = malloc(sizeof(char) * 500);
  memset(td->protein, 0, sizeof(char)*500);
  td->temp_str = malloc(sizeof(char) * 512);
  memset(td->temp_str, 0, sizeof(char)*512);

  td->insert = malloc(sizeof(int) * 100);
  td->c_delete = malloc(sizeof(int) * 100);

  int i;
  for (i=0;i<2;i++) {
    td->input_buffer[i]	= (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
    td->input_head_buffer[i] = (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
    td->output_buffer[i]	= (char **) malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
    td->aa_buffer[i]	=   (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
    td->dna_buffer[i]	= (char **)malloc(sizeof(char*) * MAX_SEQS_PER_BUFFER);
    td->acceptable_buffer[i] = (unsigned int *)malloc(sizeof(unsigned int) * MAX_SEQS_PER_BUFFER);

    int j;
    for(j=0;j<MAX_SEQS_PER_BUFFER;j++){
      td->input_buffer[i][j] = malloc(STRINGLEN);
      td->input_head_buffer[i][j] = malloc(STRINGLEN);
      td->aa_buffer[i][j] = malloc(STRINGLEN);
      td->dna_buffer[i][j] = malloc(STRINGLEN);
      td->output_buffer[i][j] = malloc(STRINGLEN);

      memset(td->output_buffer[i][j], 0, STRINGLEN);
      memset(td->input_buffer[i][j], 0, STRINGLEN);
      memset(td->dna_buffer[i][j], 0, STRINGLEN);
      memset(td->aa_buffer[i][j], 0, STRINGLEN);
      memset(td->input_head_buffer[i][j], 0, STRINGLEN );
    }
  }
}

void* writer_func(void* args) {

  int j;
  FILE* aa_outfile_fp = fopen(aa_file, "a");
  if(!aa_outfile_fp) {
    printf("ERROR: Could not open aa output file %s for writing!\n", aa_file);
    exit(0);
  } 

  while(1) {
    // write out results
    sem_wait(sema_w);

    QUEUE* temp;

    sem_wait(sema_Q);
    cutnpaste_q(&temp, DONE_Q);
    sem_post(sema_Q);

    while(temp) {

      sem_wait(sema_R);

      thread_data* td = temp->td;
      unsigned int buffer = temp->buffer;

      if (output_meta) {
        outfile_fp = fopen(out_file, "a");
        if(!outfile_fp) {
          printf("ERROR: Could not open output file %s for writing!\n", out_file);
          exit(0);
        } 
      }

      if (output_dna) {
        dna_outfile_fp = fopen(dna_file, "a");
        if (!dna_outfile_fp) {
          printf("ERROR: Could not open dna output file %s for writing!\n", dna_file);
          exit(0);
        }
      }

      num_writes++;
      for(j = 0; j < td->output_num_sequences[buffer]; j++) {
        writer_counter++;
        char *ptrc;
        if(td->aa_buffer[buffer][j][0]!=0) {
          ptrc=td->aa_buffer[buffer][j];
          while(*ptrc!='\0'){
            if(*ptrc=='\t') *ptrc='>';
            ptrc++;
          }

          fprintf(aa_outfile_fp, ">%s", td->aa_buffer[buffer][j]);
        }
        memset(td->output_buffer[buffer][j], 0, STRINGLEN);
        memset(td->aa_buffer[buffer][j], 0, STRINGLEN);
        memset(td->dna_buffer[buffer][j], 0, STRINGLEN);
      }

      if (verbose) printf("INFO: Wrote results for thread %d, buffer %d.\n", td->id, buffer );
      if (output_meta) fclose(outfile_fp);
      if (output_dna) fclose(dna_outfile_fp);

      sem_post(sema_R);
      sem_post(td->sema_w);

      temp = temp->next;
    }

    if(num_reads_flag == 1 && writer_counter ==  read_counter)   {
      sem_post(stop_sema);
      break;
    }
  }
  fclose(aa_outfile_fp);
}

void* thread_func(void *_thread_datas) {

  thread_data *td = (thread_data*)_thread_datas;
  unsigned int b = 0;
  unsigned int i;

  while(1) {
    sem_wait(td->sema_r);
    sem_wait(td->sema_w);

    sem_wait(counter_sema);
    viterbi_counter +=  td->input_num_sequences[b];
    sem_post(counter_sema);

    for (i=0; i < td->input_num_sequences[b]; i++) {
      unsigned int stringlength = strlen(td->input_buffer[b][i]);
      get_prob_from_cg(td->hmm, &train, td->input_buffer[b][i], stringlength);

      if (td->input_buffer[b][i] && td->input_head_buffer[b][i] ) {
        memset(td->aa_buffer[b][i], 0, STRINGLEN );

        viterbi(td->hmm, td->input_buffer[b][i], td->output_buffer[b][i], td->aa_buffer[b][i], td->dna_buffer[b][i], 
            td->input_head_buffer[b][i], td->wholegenome, td->format, stringlength,
            td->dna, td->dna1, td->dna_f, td->dna_f1, td->protein,
            td->insert, td->c_delete, td->temp_str);

        td->acceptable_buffer[b][i] = 1;

        sem_wait(work_sema);
        work_counter++;
        sem_post(work_sema);
      }

    }
    td->output_num_sequences[b] = td->input_num_sequences[b]; 

    if (verbose) printf("INFO: Thread %d buffer %d done work on %d sequences!\n", td->id, b, td->input_num_sequences[b]);
    sem_wait(sema_Q);
    enqueue(td, b, EMPTY_Q);
    enqueue(td, b, DONE_Q);

    sem_post(sema_Q);
    sem_post(sema_r);
    sem_post(sema_w);

    b = (b + 1) % 2;
  }
  return (void*) 0;
}
