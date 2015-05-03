#include "hmm.h"
#include "util_lib.h"
#include "fasta.h"
#include <stdio.h>
#undef READER_WRITER 
//#define READER_WRITER 

unsigned int ids; // thread id
unsigned int threadnum = 1;
unsigned int format;
unsigned int wholegenome;
unsigned int output_dna;
unsigned int output_meta;
unsigned int verbose;

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
unsigned int num_reads_flag = 0;


#ifdef MAC_SEM
sem_t *work_sema;
sem_t *stop_sema;
sem_t *counter_sema;
#else
sem_t work_sema;
sem_t stop_sema;
sem_t counter_sema;
#endif

unsigned int MAX_BYTES_PER_BUFFER;
unsigned int MAX_SEQS_PER_BUFFER;

FILE* outfile_fp;
FILE* dna_outfile_fp;

int main (int argc, char **argv) {

    int max_mem = 0;
    strncpy(train_dir, argv[0], strlen(argv[0])-4);
    strcat(train_dir, "train/");
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

    while ((c=getopt(argc, argv, "fs:m:o:w:t:p:dev")) != -1){
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
    MAX_BYTES_PER_BUFFER = max_mem*1000000/(5*2*threadnum);
    MAX_SEQS_PER_BUFFER = MAX_BYTES_PER_BUFFER/STRINGLEN;

    if (verbose)
        printf("DEBUG: Max number of sequences per thread : %d, max bytes per thread : %d\n", MAX_SEQS_PER_BUFFER, MAX_BYTES_PER_BUFFER*5);

    /* read all initial model */
    get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

    pthread_t *thread = malloc(sizeof(pthread_t*) * threadnum);
    memset(thread, 0, sizeof(pthread_t*) * threadnum);
    thread_datas = malloc(sizeof(thread_data) * threadnum);

   // memset(thread_datas, 0, sizeof(thread_data) * threadnum);

   // FILE* fp = 0;
    FASTAFILE* fp = 0;

#ifdef MAC_SEM
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

#else
    sem_init(&work_sema, 0, 1);
    sem_init(&sema_Q, 0, 1);   
    sem_init(&sema_R, 0, 1);
    sem_init(&sema_r, 0, 0);
    sem_init(&sema_w, 0, 0);
    sem_init(&stop_sema, 0, 0);
    sem_init(&counter_sema, 0, 1);
#endif

    // remove them, if they already exist
    remove(aa_file);
    if (output_meta) remove(out_file);
    if (output_dna) remove(dna_file);

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

    //fp = fopen (seq_file, "r");
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
            #ifdef MAC_SEM
                sem_post(thread_datas[j].sema_r);                    
            #else
                sem_post(&thread_datas[j].sema_r);                    
            #endif
            }
        }
    }

    if (verbose)
        printf("INFO : Initializing worker threads...\n");


    // main master loop - while we haven't exhausted reading the file yet

    while (stopped_at_fpos!=0) {
#ifdef MAC_SEM
        sem_wait(sema_r);
#else
        sem_wait(&sema_r);
#endif



#ifdef MAC_SEM
        sem_wait(sema_Q);        
        QUEUE* temp;
        cutnpaste_q(&temp, EMPTY_Q);
        sem_post(sema_Q);
#else
        sem_wait(&sema_Q);        
        QUEUE* temp;
        cutnpaste_q(&temp, EMPTY_Q);
        sem_post(&sema_Q);
#endif

        while(temp) {
#ifdef MAC_SEM
            sem_wait(sema_R);
#else
            sem_wait(&sema_R);
#endif

            //stopped_at_fpos = read_file_into_buffer(fp, stopped_at_fpos, temp->td, temp->buffer);
            stopped_at_fpos = read_seq_into_buffer(fp,  temp->td, temp->buffer);

#ifdef MAC_SEM
            sem_post(sema_R);
#else
            sem_post(&sema_R);
#endif


#ifdef VERBOSE
            printf("READER subsequent reads into %d seqs into id %d\n", temp->td->num_sequences[temp->buffer], temp->td->id);
#endif 

#ifdef MAC_SEM
            sem_post(temp->td->sema_r);
#else
            sem_post(&temp->td->sema_r);
#endif
            temp = temp->next;
        }
    }
    CloseFASTA(fp);

    if (verbose)
        printf("INFO : Finished handing out all the work...\n");

    //total_reads = read_counter;
   // printf("Done reading %ld \n", total_reads);

    num_reads_flag =1;

#ifdef MAC_SEM
    sem_wait(stop_sema);
#else
    sem_wait(&stop_sema);
#endif

#ifdef MAC_SEM
    sem_unlink("/work_sema");
    sem_unlink("/sema_Q");
    sem_unlink("/sema_R");
    sem_unlink("/sema_r");
    sem_unlink("/sema_w");
    sem_unlink("/stop_sema");
    sem_unlink("/COUNTER_SEMA");

    char name[40];
    for(j=0; j<threadnum; j++) {
       sprintf(name, "/sema_r%d", j);
       sem_unlink(name);

       sprintf(name, "/sema_w%d", j);
       sem_unlink(name);
    }

#endif

    printf("Run finished with %d threads.\n", threadnum);
    printf("#Seq read :%ld, %ld  #seq viterbied : %ld  # seq written :  %ld \n", read_counter, read_counter1,  viterbi_counter, writer_counter);
}
char mystring[STRINGLEN];
char complete_sequence[STRINGLEN];


int read_seq_into_buffer(FASTAFILE* ffp, thread_data* thread_data, unsigned int buf) {

  char *seq;
  char *name;
  int   L;
  int status ;
  int count=0;

  while ( (count < MAX_SEQS_PER_BUFFER) && (status = ReadFASTA(ffp, &seq, &name, &L)) ==1 ) {
      //printf(">%s\n", name);
      //printf("%s\n",  seq);
      strcpy(thread_data->input_head_buffer[buf][count], name);
      strcpy(thread_data->input_buffer[buf][count], seq);
      read_counter++;
    //  free(seq);
     // free(name);
      count++;
  }

  thread_data->input_num_sequences[buf] = count;
  read_counter1 += thread_data->input_num_sequences[buf];

#ifdef READER_WRITER
    printf("READER : %d - %d  %d\n", thread_data->id, buf, count);
#endif


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

   // char* mystring = (char*) malloc(sizeof(char) * STRINGLEN);
   // char* complete_sequence = (char*) malloc(sizeof(char) * STRINGLEN);
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
                //if (strlen(complete_sequence) > 70) { read_counter++;}
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

#ifdef READER_WRITER
    printf("READER : %d - %d  %ld\n", thread_data->id, buf, count);
#endif


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

#ifdef MAC_SEM
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


   // sem_post(td->sema_w);
  //  sem_post(td->sema_w);

#else
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

        /*   memset(td->input_buffer[i], 0, sizeof(char*) * MAX_SEQS_PER_BUFFER);
             memset(td->input_head_buffer[i], 0, sizeof(char*) * MAX_SEQS_PER_BUFFER);
             memset(td->output_buffer[i], 0, sizeof(char*) * MAX_SEQS_PER_BUFFER);
             memset(td->aa_buffer[i], 0, sizeof(char*) * MAX_SEQS_PER_BUFFER);
             memset(td->dna_buffer[i], 0, sizeof(char*) * MAX_SEQS_PER_BUFFER);
             memset(td->acceptable_buffer[i], 0, sizeof(unsigned int) * MAX_SEQS_PER_BUFFER);
             */

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
#ifdef MAC_SEM
        sem_wait(sema_w);
#else
        sem_wait(&sema_w);
#endif

        QUEUE* temp;

#ifdef MAC_SEM
        sem_wait(sema_Q);
#else
        sem_wait(&sema_Q);
#endif
        cutnpaste_q(&temp, DONE_Q);
#ifdef MAC_SEM
        sem_post(sema_Q);
#else
        sem_post(&sema_Q);
#endif

        while(temp) {

#ifdef MAC_SEM
            sem_wait(sema_R);
#else
            sem_wait(&sema_R);
#endif

#ifdef VERBOSE
            printf("WRITER : tid %d temp %p buffer %d\n", temp->td->id, temp, temp->buffer); 
#endif
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
                //fprintf(aa_outfile_fp, ">xx\n");
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
                else {
                  //fprintf(aa_outfile_fp, ">%s", td->aa_buffer[buffer][j]);
                }
                if (td->acceptable_buffer[j]) {

                   // fprintf(aa_outfile_fp, "%s", td->aa_buffer[buffer][j]);
                    //  if (output_meta) fprintf(outfile_fp, "%s", td->output_buffer[j]);
                    // if (output_dna) fprintf(dna_outfile_fp, "%s", td->dna_buffer[j]);
                }
                memset(td->output_buffer[buffer][j], 0, STRINGLEN);
                memset(td->aa_buffer[buffer][j], 0, STRINGLEN);
                memset(td->dna_buffer[buffer][j], 0, STRINGLEN);
            }



            if (verbose) printf("INFO: Wrote results for thread %d, buffer %d.\n", td->id, buffer );

            if (output_meta) fclose(outfile_fp);
            if (output_dna) fclose(dna_outfile_fp);

            //td->num_sequences[buffer] = 0;

#ifdef MAC_SEM
            sem_post(sema_R);
            sem_post(td->sema_w);
#else
            sem_post(&sema_R);
            sem_post(&td->sema_w);
#endif

            temp = temp->next;
            /*
               old_temp = temp;
               temp = temp->next ? temp->next : 0;
               free(old_temp); 
               */
        }
    //    printf("Total reads = %ld   total writes = %ld\n", num_reads, num_writes);
#ifdef READER_WRITER
       printf("Total reads = %ld   total writes = %ld\n", read_counter, writer_counter);
#endif
     
        if(num_reads_flag == 1 && writer_counter ==  read_counter)   {

#ifdef MAC_SEM
          printf("TERMINATING\n");
          sem_post(stop_sema);
#else
          sem_post(&stop_sema);
#endif
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
#ifdef MAC_SEM
        sem_wait(td->sema_r);
        sem_wait(td->sema_w);
#else
        sem_wait(&td->sema_r);
        sem_wait(&td->sema_w);
#endif


#ifdef READER_WRITER
        printf("   WORKER : %d - %d   WORK : %d\n", td->id, b,  td->input_num_sequences[b]);
#endif



#ifdef MAC_SEM
        sem_wait(counter_sema);
        viterbi_counter +=  td->input_num_sequences[b];
        sem_post(counter_sema);
#else
        sem_wait(&counter_sema);
        viterbi_counter +=  td->input_num_sequences[b];
        sem_post(&counter_sema);
#endif

        //viterbi_counter +=  td->num_sequences[b];
        for (i=0; i < td->input_num_sequences[b]; i++) {
            unsigned int stringlength = strlen(td->input_buffer[b][i]);
            get_prob_from_cg(td->hmm, &train, td->input_buffer[b][i], stringlength);
   
            if(td->input_buffer[b][i] == 0 || td->input_head_buffer[b][i] == 0 ) {
              printf("%s\n",td->input_buffer[b][i]);  
              printf("%s\n",td->input_head_buffer[b][i]);
            }


            if (td->input_buffer[b][i] != 0 && td->input_head_buffer[b][i] != 0 ) {
                memset(td->aa_buffer[b][i], 0, STRINGLEN );
                //memset(td->dna_buffer[b][i], 0, STRINGLEN );


                viterbi(td->hmm, td->input_buffer[b][i], td->output_buffer[b][i], td->aa_buffer[b][i], td->dna_buffer[b][i], 
                        td->input_head_buffer[b][i], td->wholegenome, td->format, stringlength,
                        td->dna, td->dna1, td->dna_f, td->dna_f1, td->protein,
                        td->insert, td->c_delete, td->temp_str);

                td->acceptable_buffer[b][i] = 1;
               // strcpy(td->output_buffer[b][i], td->aa_buffer[b][i]);
                //printf("%s",td->output_buffer[b][i]);

#ifdef MAC_SEM
                sem_wait(work_sema);
                work_counter++;
                sem_post(work_sema);
#else
                sem_wait(&work_sema);
                work_counter++;
                sem_post(&work_sema);
#endif
            }

        }
        td->output_num_sequences[b] = td->input_num_sequences[b]; 

        if (verbose) printf("INFO: Thread %d buffer %d done work on %d sequences!\n", td->id, b, td->input_num_sequences[b]);
#ifdef MAC_SEM
        sem_wait(sema_Q);
#else
        sem_wait(&sema_Q);
#endif
        enqueue(td, b, EMPTY_Q);
        enqueue(td, b, DONE_Q);

#ifdef MAC_SEM
        sem_post(sema_Q);
        sem_post(sema_r);
        sem_post(sema_w);
#else
        sem_post(&sema_Q);
        sem_post(&sema_r);
        sem_post(&sema_w);
#endif


        b = (b + 1) % 2;
    }
    return (void*) 0;
}
