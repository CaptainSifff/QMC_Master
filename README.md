This is the master process of my Monte-Carlo code.
It waits for connections from clients and givs them the parameters for the simulation.
After every measurement the data is collected by the server. This way we can resist the loss of clients
due to some external circumstances.
The error analysis and storage of data is only performed by the master.