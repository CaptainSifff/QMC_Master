My Monte Carlo Code is built as a client-server architecture.
This is the master process of my Monte Carlo code.
It waits for connections from clients and distributes to them the parameters for the simulation.
After every measurement the data is collected by the server. This way we are able to withstand the loss of clients
due to some external circumstances without interrupting the work of other nodes.
The error analysis and storage of data is only performed by the master.

The CMakelists expects in the mpi Subfolder a symlink onto the cpp files in the SDL subfolder.
