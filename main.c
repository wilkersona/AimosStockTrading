#include "clockcycle.h"
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define TOTAL_BUYERS 100
#define TOTAL_STOCKS 10
#define NUM_T 2

struct Buyer {
    int id;
    int portfolio[TOTAL_STOCKS];
    int cash;
    int buy_strat;
    int sell_strat;
    int commitment;
};

struct Stock {
    int id;
    float value_m2;
    float value_m1;
    float value;
    float dividend;
    float volatility;
};

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
    for (int i=0; i<buyers_per; i++) {
        // PUT LOGIC TO RANDOMLY INITIALIZE BUYER HERE
        // CANT BE SEEDED OTHERWISE ALL RANKS WILL HAVE SAME BUYERS
        struct Buyer new_buyer;
        new_buyer.id = (rank-1)*buyers_per + i;
        new_buyer.cash = 100;
        new_buyer.buy_strat = 0;
        new_buyer.sell_strat = 0;
        new_buyer.commitment = 0;
        for (int j=0; j<TOTAL_STOCKS; j++) {
            new_buyer.portfolio[j] = (j == (rank-1)*buyers_per + i) ? 1 : 0;
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
        // PUT LOGIC TO RANDOMLY INITIALIZE STOCK HERE
        struct Stock new_stock;
        new_stock.id = i;
        new_stock.value_m2 = 10.0;
        new_stock.value_m1 = 10.0;
        new_stock.value = i+1;
        new_stock.dividend = 0.01;
        new_stock.volatility = 0.05;
        // Make sure to leave a buffer of <stocks_per> at the front of the list - this is the section the I/O rank "contributes"
        out[i+stocks_per] = new_stock;
    }
    // SORT STOCKS BY DIVIDEND
    return out;
}

void save_stocks(struct Stock* stocks, int stocks_per) {
    /*
    * I/O call that saves the stock data
    *
    * Arguments-
    * stocks:               list of all stocks, padded with <stocks_per> blank entries for the I/O rank
    * stocks_per:           number of stocks per rank
    */
    printf("Stock values:\n");
    for (int i=stocks_per; i<TOTAL_STOCKS+stocks_per; i++) {
        // PUT LOGIC TO SAVE STOCK VALUE HERE
        printf("%i, %.1f ", stocks[i].id, stocks[i].value);
    }
    printf("\n");
}

void save_buyers(float* values, int buyers_per) {
    /*
    * I/O call that saves the buyer data
    *
    * Arguments-
    * values:               list of paired values for all buyers id0, value0, id1, value1, etc...
    * buyers_per:           number of buyers per rank
    */
    printf("Buyer values:\n");
    for (int i=buyers_per; i<TOTAL_BUYERS+buyers_per; i++) {
        // PUT LOGIC TO SAVE STOCK VALUE HERE
        // NOTICE EACH ID IS IN i*2 AND EACH VALUE IS IN i*2+1
        printf("%i, %.1f ", (int) values[i*2], values[i*2+1]);
    }
    printf("\n");
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
        // Might eventually include dividend as well
        out = out + buyer.portfolio[i]*stocks[i+stocks_per].value;
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
    * (float*):             list of pairs of ids and values. For buyer i the id is in index i*2 and the value in i*2+1
    */
    for (int i=0; i<buyers_per; i++) {
        // Assign the id and the total value of each buyer
        buyer_values[i*2] = (float) buyers[i].id;
        buyer_values[i*2+1] = get_total_value(buyers[i], stocks, stocks_per);
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
        for (int j=0; j<TOTAL_STOCKS; j++) {
            // Chunk is the rank that receives the send
            chunk = j / stocks_per + 1;
            // Index is the position of the stock in the chunk
            index = j % stocks_per;

            // PUT LOGIC HERE TO CHECK STOCK DATA AGAINST STRATEGY
            // IF STOCK PASSES CREATE SEND
            // UPDATE BUYER CASH AND PORTFOLIO
            // STOCK J IS IN buyer.portfolio[j] AND stocks[j+stocks_per]
            if (i == 0 && j == 0) {
                buyers[i].portfolio[j] = buyers[i].portfolio[j] + 2;
                action[0] = 1;
                action[1] = 2;
                action[2] = index;
                MPI_Isend(action, 3, MPI_INT, chunk, 0, MPI_COMM_WORLD, &request);
            } else if (i == 1 && j == 1) {
                buyers[i].portfolio[j] = buyers[i].portfolio[j] - 1;
                action[0] = 0;
                action[1] = 1;
                action[2] = index;
                MPI_Isend(action, 3, MPI_INT, chunk, 0, MPI_COMM_WORLD, &request);
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
    // Probe for first message
    MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
    while (flag == 1) {
        // While there are messages, receive them
        MPI_Recv(action, 3, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        // PUT LOGIC HERE TO COMPLETE ACTION BASED ON EVENT
        if (action[0] == 0) { // Sell
            stocks[my_chunk + action[2]].value--;
        } else { // Buy
            stocks[my_chunk + action[2]].value++;
        }
        // Probe for next message
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
    }
}

int main(int argc, char *argv[]) {

    // Initialize MPI
    int npes, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

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
        buyer_values = (float*) calloc(buyers_per*2, sizeof(float));
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
        //buyer_values = calloc(buyers_per*2, sizeof(float));
        all_values = calloc(TOTAL_BUYERS*2+buyers_per*2, sizeof(float));
    } else {
        // Worker ranks generate their buyers' values
        buyer_values = value_chunk(buyer_values, buyers, all_stocks, buyers_per, stocks_per);
    }

    // Gather initial buyer values to the I/O rank for saving
    MPI_Gather(
        myrank == 0 ? MPI_IN_PLACE : buyer_values,
        buyers_per*2,
        MPI_FLOAT,
        all_values,
        buyers_per*2,
        MPI_FLOAT,
        0,
        MPI_COMM_WORLD
    );

    if (myrank == 0) {
        // Save initial buyer values and intial stock values
        save_buyers(all_values, buyers_per);
        save_stocks(all_stocks, stocks_per);
    }

    // Determines the starting index of a workers stock chunk, INCLUDES PADDED VALUES FOR RANK 0
    int my_chunk = myrank*stocks_per;

    // Run simulation
    int t = 0;
    while (t < NUM_T) {
        
        if (myrank != 0) {
            // All worker ranks update their buyers and generate their values
            update_buyers(buyers, all_stocks, buyers_per, stocks_per);
            buyer_values = value_chunk(buyer_values, buyers, all_stocks, buyers_per, stocks_per);
        }

        // Gathering buyer values serves as a barrier for all ranks
        MPI_Gather(
            myrank == 0 ? MPI_IN_PLACE : buyer_values,
            buyers_per*2,
            MPI_FLOAT,
            all_values,
            buyers_per*2,
            MPI_FLOAT,
            0,
            MPI_COMM_WORLD
        );

        if (myrank == 0) {
            // I/O rank saves buyer values while
            save_buyers(all_values, buyers_per);
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
            save_stocks(all_stocks, stocks_per);
        }

        t++;
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