#include "clockcycle.h"
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TOTAL_BUYERS 100
#define TOTAL_STOCKS 10
#define NUM_T 2
#define RAND_SEED 42

struct Buyer {
    int id;
    int portfolio[TOTAL_STOCKS];
    int cash;
    int allowance;
    int buy_strat;
    int sell_strat;
    float commitment;
};

struct Stock {
    int id;
    float value_m2;
    float value_m1;
    float value;
    float dividend;
    float volatility;
};

int stock_comparator(const void *v1, const void *v2)
{
    const struct Stock *p1 = (struct Stock *) v1;
    const struct Stock *p2 = (struct Stock *) v2;
    if (p1->dividend > p2->dividend)
        return -1;
    else if (p1->dividend < p2->dividend)
        return +1;
    else
        return 0;
}

struct Buyer* init_buyers(int rank, int buyers_per) {
    /*
    * Initializes the buyers for the given rank
    * 
    * Arguments-
    * rank:                 the current rank
    * buyers_per:           the number of buyers per worker rank
    * 
    * Returns-
    * (struct Buyer*):      list of buyers randomly initialized
    */
    struct Buyer* out = (struct Buyer*) calloc(buyers_per, sizeof(struct Buyer)); // Allocate space for list
    int non_overlap[3];
    //Init randomness
    srand(RAND_SEED + rank);
    
    for (int i=0; i<buyers_per; i++) {
    	//Logic to create new buyer
        struct Buyer new_buyer;
        new_buyer.id = (rank-1)*buyers_per + i;
        new_buyer.cash = 100;
        new_buyer.allowance = 5;
        new_buyer.buy_strat = rand()%4;		//KEY: 0: Rise, 1: Fall, 2: Peak, 3: Dip
        for (int strat=0; strat<3; strat++) {
            if (strat >= new_buyer.buy_strat) {
                non_overlap[strat] = strat + 1;
            } else {
                non_overlap[strat] = strat;
            }
        }
        new_buyer.sell_strat = non_overlap[rand()%3];
        new_buyer.commitment = (1 + rand()%100)/100.0; //Instead of half, full, etc. we can mark this as a percent?
        for (int j=0; j<TOTAL_STOCKS; j++) {
            //new_buyer.portfolio[j] = (rand()%5)/4; //20% chance each stock is in portfolio
            //new_buyer.portfolio[j] = (j == (rank-1)*buyers_per + i) ? 1 : 0;
            new_buyer.portfolio[j] = 0;
        }
        out[i] = new_buyer;
    }
    return out;
}

struct Stock* init_stocks(int stocks_per) {
    /*
    * Initializes all the stocks for the simulation
    *
    * Arguments-
    * stocks_per:           number of stocks per rank, used to pad the list for the I/O rank
    *
    * Returns-
    * (struct Stock*):      list of stocks randomly initialized, SORTED BY DIVIDEND
    */
    struct Stock* out = (struct Stock*) calloc(TOTAL_STOCKS + stocks_per, sizeof(struct Stock)); // Allocate space for list
    for (int i=0; i<TOTAL_STOCKS; i++) {
        //Init Randomness
    	srand(RAND_SEED + i);

    	//Stock variables initialization
        struct Stock new_stock;
        new_stock.id = i;
        new_stock.value_m2 = (50 + rand()%100)/10.0;
        new_stock.value_m1 = (50 + rand()%100)/10.0;
        new_stock.value = (50 + rand()%100)/10.0;
        new_stock.dividend = (5 + rand()%10)/1000.0;
        new_stock.volatility = (1 + rand()%10)/100.0;
        // Make sure to leave a buffer of <stocks_per> at the front of the list - this is the section the I/O rank "contributes"
        out[i+stocks_per] = new_stock;
    }
    // SORT STOCKS BY DIVIDEND
    qsort(out+stocks_per, TOTAL_STOCKS, sizeof(struct Stock), stock_comparator);
    return out;
}

unsigned int save_stocks(struct Stock* stocks, int stocks_per, MPI_File file) {
    /*
    * I/O call that saves the stock data
    *
    * Arguments-
    * stocks:               list of all stocks, padded with <stocks_per> blank entries for the I/O rank
    * stocks_per:           number of stocks per rank
    *
    * Returns-
    * (unsigned int):       total number of bytes written
    */
    MPI_Status status;
    MPI_Offset offset;
    MPI_File_open(MPI_COMM_SELF, "stocks", MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_File_get_position(file, &offset);
    MPI_File_seek(file, offset, MPI_SEEK_SET);
    float* values = (float*) calloc(TOTAL_STOCKS, sizeof(float));
    for (int i=0; i<TOTAL_STOCKS; i++) {
        values[i] = stocks[i+stocks_per].value;
    }
    MPI_File_write(file, values, TOTAL_STOCKS, MPI_FLOAT, &status);
    MPI_File_close(&file);
    free(values);
    return TOTAL_STOCKS*sizeof(float);
}

unsigned int save_buyers(float* values, int buyers_per, MPI_File file) {
    /*
    * I/O call that saves the buyer data
    *
    * Arguments-
    * values:               list of values for all buyers
    * buyers_per:           number of buyers per rank
    *
    * Returns-
    * (unsigned int):       total number of bytes written
    */
    MPI_Status status;
    MPI_Offset offset;
    MPI_File_open(MPI_COMM_SELF, "buyers", MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_File_get_position(file, &offset);
    MPI_File_seek(file, offset, MPI_SEEK_SET);
    MPI_File_write(file, values+buyers_per, TOTAL_BUYERS, MPI_FLOAT, &status);
    MPI_File_close(&file);
    return TOTAL_BUYERS*sizeof(float);
}

void save_strategies(float* strategies, int buyers_per) {
    /*
    * I/O call that saves the buyer strategies, only called once
    *
    * Arguments-
    * strategies:           list of values for all buyers
    * buyers_per:           number of buyers per rank
    */
    MPI_Status status;
    MPI_File strat_file;
    MPI_File_open(MPI_COMM_SELF, "strats", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &strat_file);
    MPI_File_write(strat_file, strategies+buyers_per*3, TOTAL_BUYERS*3, MPI_FLOAT, &status);
    MPI_File_close(&strat_file);
}

float get_total_value(struct Buyer buyer, struct Stock* stocks, int stocks_per) {
    /*
    * Gets the total value of a buyer's portfolio and cash
    * 
    * Arguments-
    * buyer:                buyer you want the value of
    * stocks:               list of all stocks, padded with <stocks_per> blank entries for the I/O rank
    * stocks_per:           number of stocks per rank
    *
    * Returns-
    * (float):              sum of the buyer's cash with the value of the buyer's portfolio
    */
    float out = buyer.cash;
    for (int i=0; i<TOTAL_STOCKS; i++) {
        // Value of the portfolio if the buyer sold all stocks at their current value
        out += buyer.portfolio[i]*stocks[i+stocks_per].value + buyer.portfolio[i]*stocks[i+stocks_per].value*stocks[i+stocks_per].dividend;
    }
    return out;
}

float* value_chunk(float* buyer_values, struct Buyer* buyers, struct Stock* stocks, int buyers_per, int stocks_per) {
    /*
    * Gets the list of id, value pairs for a list of buyers
    *
    * Arguments-
    * buyer_values:         allocated list of memory to place data in
    * buyers:               list of the current rank's buyers
    * stocks:               list of all stocks, padded with <stocks_per> blank entries for the I/O rank
    * buyers_per:           number of buyers per rank
    * stocks_per:           number of stocks per rank
    *
    * Returns-
    * (float*):             list values
    */
    for (int i=0; i<buyers_per; i++) {
        // Assign the id and the total value of each buyer
        buyer_values[i] = get_total_value(buyers[i], stocks, stocks_per);
    }
    return buyer_values;
}

void update_buyers(struct Buyer* buyers, struct Stock* stocks, int buyers_per, int stocks_per) {
    /*
    * For a chunk of buyers determines and sends off the stock events that the buyer completes
    * Updates the values of the buyer according to their actions
    *
    * Arguments-
    * buyers:               list of the current rank's buyers
    * stocks:               list of all stocks, padded with <stocks_per> blank entries for the I/O rank
    * buyers_per:           number of buyers per rank
    * stocks_per:           number of stocks per rank
    */
    MPI_Request request;
    int chunk;
    int index;
    int action[3];      // action[0] = {0, 1} for sell/buy, action[1] = # bought/sold, action[2] = index of stock in chunk
    for (int i=0; i<buyers_per; i++) {
        // Allowance
        buyers[i].cash += buyers[i].allowance;

        // SELLING
        for (int j=0; j<TOTAL_STOCKS; j++) {
            // Chunk is the rank that receives the send
            chunk = j / stocks_per + 1;
            // Index is the position of the stock in the chunk
            index = j % stocks_per;

            // Selling based on sell strategy
            if((buyers[i].sell_strat == 0 && stocks[j+stocks_per].value >= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 >= stocks[j+stocks_per].value_m2)
            || (buyers[i].sell_strat == 1 && stocks[j+stocks_per].value <= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 <= stocks[j+stocks_per].value_m2)
            || (buyers[i].sell_strat == 2 && stocks[j+stocks_per].value <= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 >= stocks[j+stocks_per].value_m2)
            || (buyers[i].sell_strat == 3 && stocks[j+stocks_per].value >= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 <= stocks[j+stocks_per].value_m2)) { 
                
                action[0] = 0;
                action[1] = (int) ceil((buyers[i].portfolio[j] * buyers[i].commitment)); 
                buyers[i].cash += stocks[j+stocks_per].value*action[1];
                buyers[i].portfolio[j] -= action[1];

                action[2] = index;
                // If any stocks are sold, send the action
                if (action[1] != 0)
                    MPI_Isend(action, 3, MPI_INT, chunk, 0, MPI_COMM_WORLD, &request);
            }
        }

        // BUYING
        for (int j=0; j<TOTAL_STOCKS; j++) {
            // Chunk is the rank that receives the send
            chunk = j / stocks_per + 1;
            // Index is the position of the stock in the chunk
            index = j % stocks_per;

            // STOCK J IS IN buyer.portfolio[j] AND stocks[j+stocks_per]

            // Buying based on buy strategy
            if((buyers[i].buy_strat == 0 && stocks[j+stocks_per].value >= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 >= stocks[j+stocks_per].value_m2)
            || (buyers[i].buy_strat == 1 && stocks[j+stocks_per].value <= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 <= stocks[j+stocks_per].value_m2)
            || (buyers[i].buy_strat == 2 && stocks[j+stocks_per].value <= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 >= stocks[j+stocks_per].value_m2)
            || (buyers[i].buy_strat == 3 && stocks[j+stocks_per].value >= stocks[j+stocks_per].value_m1 && stocks[j+stocks_per].value_m1 <= stocks[j+stocks_per].value_m2)) { 
                
                action[0] = 1;
                action[1] = 0;
                float original_cash = buyers[i].cash;
                while (buyers[i].cash >= 0) {
                    if (buyers[i].cash - stocks[j+stocks_per].value >= original_cash*(1.0-buyers[i].commitment)) {
                        buyers[i].cash -= stocks[j+stocks_per].value;
                        buyers[i].portfolio[j]++;
                        action[1]++;
                    }
                    else
                        break;
                }
                
                action[2] = index;
                // If any stocks are bought, send the action
                if (action[1] != 0)
                    MPI_Isend(action, 3, MPI_INT, chunk, 0, MPI_COMM_WORLD, &request);
                break;
            }
        }
    }
}

void update_stocks(int my_chunk, struct Stock* stocks, int stocks_per) {
    /*
    * For a chunk of stocks pull events corresponding to the chunk and complete the action
    * 
    * Arguments-
    * my_chunk:             index at which the current rank's chunk of stocks begins
    * stocks:               list of all stocks, padded with <stocks_per> blank entries for the I/O rank
    * stocks_per:           number of stocks per rank
    */
    int action[3];
    int flag = 0;
    MPI_Status status;
    // Save value_m2 and value_m1
    for (int i=0; i<stocks_per; i++) {
        stocks[my_chunk + i].value_m2 = stocks[my_chunk + i].value_m1;
        stocks[my_chunk + i].value_m1 = stocks[my_chunk + i].value;         // USE THIS TO SCALE VOLATILITY
    }
    // Probe for first message
    MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
    while (flag == 1) {
        // While there are messages, receive them
        MPI_Recv(action, 3, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        // PUT LOGIC HERE TO COMPLETE ACTION BASED ON EVENT
        if (action[0] == 0) { // Sell
            stocks[my_chunk + action[2]].value -= stocks[my_chunk + action[2]].value_m1*stocks[my_chunk + action[2]].volatility*action[1];
            if (stocks[my_chunk + action[2]].value <= 1) {
                stocks[my_chunk + action[2]].value = 1;
            }
        } else { // Buy
            stocks[my_chunk + action[2]].value += stocks[my_chunk + action[2]].value_m1*stocks[my_chunk + action[2]].volatility*action[1];
        }
        // Probe for next message
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
    }
}

int main(int argc, char *argv[]) {

    // Initialize clock
    unsigned long long start_time;
    unsigned long long end_time;
    unsigned long long start_time_IO;
    unsigned long long end_time_IO;
    unsigned long long start_time_perceived;
    unsigned long long end_time_perceived;

    // Initialize MPI
    int npes, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Define files
    MPI_File buyer_file;
    MPI_File stock_file;
    if (myrank == 0) {
        MPI_File_open(MPI_COMM_SELF, "buyers", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &buyer_file);
        MPI_File_open(MPI_COMM_SELF, "stocks", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &stock_file);
        MPI_File_close(&buyer_file);
        MPI_File_close(&stock_file);
    }

    // Define new MPI datatype for sending Stock structs
    const int n_members = 6;
    int blocklengths[6] = {1, 1, 1, 1, 1, 1};
    MPI_Datatype types[6] = {MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_stock_type;
    MPI_Aint offsets[6];

    offsets[0] = offsetof(struct Stock, id);
    offsets[1] = offsetof(struct Stock, value_m2);
    offsets[2] = offsetof(struct Stock, value_m1);
    offsets[3] = offsetof(struct Stock, value);
    offsets[4] = offsetof(struct Stock, dividend);
    offsets[5] = offsetof(struct Stock, volatility);

    MPI_Type_create_struct(n_members, blocklengths, offsets, types, &mpi_stock_type);
    MPI_Type_commit(&mpi_stock_type);

    // Remove one process for I/O
    npes--;

    // Split up buyers and stocks, TOTAL_BUYERS and TOTAL_STOCKS should be divisible by (total number of ranks-1)
    int buyers_per = TOTAL_BUYERS/npes;
    int stocks_per = TOTAL_STOCKS/npes;

    if (buyers_per*npes != TOTAL_BUYERS || stocks_per*npes != TOTAL_STOCKS) {
        printf("Number of buyers or stocks not evenly divisible by processes: %i\n", npes);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Initialize buyers, each worker rank initializes their buyers separately
    struct Buyer* buyers;
    if (myrank != 0) {
        buyers = init_buyers(myrank, buyers_per);
    }
    // Initialize local buyer value chunk and I/O all values
    float* buyer_values;
    float* all_values;
    float* all_strats;
    float* local_strats;
    // Initialize local stock chunk and all stocks list
    struct Stock* stock_chunk;
    struct Stock* all_stocks;

    if (myrank == 0) {
        // I/O rank generates the stocks
        all_stocks = init_stocks(stocks_per);
    } else {
        // Worker ranks allocate memory to place the stocks and buyer values in
        stock_chunk = (struct Stock*) calloc(TOTAL_STOCKS, sizeof(struct Stock));
        all_stocks = (struct Stock*) calloc(TOTAL_STOCKS+stocks_per, sizeof(struct Stock));
        buyer_values = (float*) calloc(buyers_per, sizeof(float));
        local_strats = (float*) calloc(buyers_per*3, sizeof(float));
    }

    // Initial stock values need to be broadcasted to all ranks
    MPI_Bcast(
        all_stocks,
        TOTAL_STOCKS+stocks_per,
        mpi_stock_type,
        0,
        MPI_COMM_WORLD
    );

    if (myrank == 0) {
        // I/O rank allocates memory for all buyer values 
        all_values = calloc(TOTAL_BUYERS+buyers_per, sizeof(float));
        all_strats = calloc((TOTAL_BUYERS+buyers_per)*3, sizeof(float));
    } else {
        // Worker ranks generate their buyers' values
        buyer_values = value_chunk(buyer_values, buyers, all_stocks, buyers_per, stocks_per);
        for (int i=0; i<buyers_per; i++) {
            local_strats[i*3] = buyers[i].buy_strat;
            local_strats[i*3+1] = buyers[i].sell_strat;
            local_strats[i*3+2] = buyers[i].commitment;
        }
    }

    // Gather initial buyer values to the I/O rank for saving
    MPI_Gather(
        myrank == 0 ? MPI_IN_PLACE : buyer_values,
        buyers_per,
        MPI_FLOAT,
        all_values,
        buyers_per,
        MPI_FLOAT,
        0,
        MPI_COMM_WORLD
    );

    // Gather buyer strategies to the I/O rank
    MPI_Gather(
        myrank == 0 ? MPI_IN_PLACE : local_strats,
        buyers_per*3,
        MPI_FLOAT,
        all_strats,
        buyers_per*3,
        MPI_FLOAT,
        0,
        MPI_COMM_WORLD
    );

    if (myrank == 0) {
        // Save initial buyer values and intial stock values
        save_buyers(all_values, buyers_per, buyer_file);
        save_stocks(all_stocks, stocks_per, stock_file);
        save_strategies(all_strats, buyers_per);
        free(all_strats);
    } else {
        free(local_strats);
    }

    // Determines the starting index of a workers stock chunk, INCLUDES PADDED VALUES FOR RANK 0
    int my_chunk = myrank*stocks_per;

    // Run simulation
    int t = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = clock_now();
    float total_IO_time = 0;
    float perceived_IO_time = 0;
    unsigned int IO_size = 0;
    while (t < NUM_T) {
        
        
        if (myrank != 0) {
            // All worker ranks update their buyers and generate their values
            update_buyers(buyers, all_stocks, buyers_per, stocks_per);
            buyer_values = value_chunk(buyer_values, buyers, all_stocks, buyers_per, stocks_per);
        }
        
        // Gathering buyer values serves as a barrier for all ranks
        start_time_perceived=clock_now();
        MPI_Gather(
            myrank == 0 ? MPI_IN_PLACE : buyer_values,
            buyers_per,
            MPI_FLOAT,
            all_values,
            buyers_per,
            MPI_FLOAT,
            0,
            MPI_COMM_WORLD
        );
        end_time_perceived=clock_now();
        perceived_IO_time += ((float)(end_time_perceived - start_time_perceived)) / 512000000.0;

        if (myrank == 0) {
            // I/O rank saves buyer values while
            start_time_IO=clock_now();
            IO_size += save_buyers(all_values, buyers_per, buyer_file);
            end_time_IO=clock_now();
            total_IO_time += ((float)(end_time_IO - start_time_IO)) / 512000000.0;
        } else {
            // worker ranks update their stock chunk
            update_stocks(my_chunk, all_stocks, stocks_per);
        }
        
        if (myrank != 0) {
            // Copy updated stock values into stock_chunk
            for (int i=0; i<stocks_per; i++) {
                stock_chunk[i] = all_stocks[i+my_chunk];
            }
        }

        // Gather updated stock chunks into fully updated "global" stock array
        // Serves as a barrier for all ranks (so buyer I/O will be done and all stock chunks will be updated)
        MPI_Allgather(
            myrank == 0 ? MPI_IN_PLACE : stock_chunk,
            stocks_per,
            mpi_stock_type,
            all_stocks,
            stocks_per,
            mpi_stock_type,
            MPI_COMM_WORLD
        );

        if (myrank == 0) {
            // I/O rank saves stock data, worker ranks can continue in the loop and start updating their buyers
            start_time_IO=clock_now();
            IO_size += save_stocks(all_stocks, stocks_per, stock_file);
            end_time_IO=clock_now();
            total_IO_time += ((float)(end_time_IO - start_time_IO)) / 512000000.0; 
        }

        t++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = clock_now();
    float time_in_secs_overall = ((float)(end_time - start_time)) / 512000000.0;
    float bandwidth = ((float) IO_size) / total_IO_time;
    if (myrank == 0) {
        printf("Overall computation loop took %f seconds\n", time_in_secs_overall);
        printf("I/O took %f seconds\n", total_IO_time);
        printf("Sent %u bytes at a bandwidth of %f bytes per second\n", IO_size, bandwidth);
    } else if (myrank == 1) {
        printf("Perceived IO took %f seconds\n", perceived_IO_time);
    }

    // Frees
    if (myrank != 0) {
        free(buyers);
        free(stock_chunk);
        free(buyer_values);
    } else {
        free(all_values);
    }
    free(all_stocks);
    MPI_Type_free(&mpi_stock_type);
    MPI_Finalize();
    return 0;
}