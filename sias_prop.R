library(igraph)
library(tidyverse)

# we will generate the data using the propogation process described in the paper

# first generate a graph
set.seed(42)
n <- 100
p <- 0.1 # how connected the nodes are
G <- erdos.renyi.game(n,p,"gnp")
G_M <- as_adjacency_matrix(G)

# randomly pick source
s <- sample(1:n,1)

# pick sim variables
t_limit <- 100 # time_limit
alpha_ij <- matrix(runif(n*n,0,20),n,n) # transmission parameter
beta_i <- runif(n,0,100) # incubation parameter

# set timer
t <- 0

# create dataset for keeping up with sim
nodes <- tibble(node=1:n, # node id
                infectious=rep(F,n), # whether or not infectious
                incubating_until=rep(-Inf,n), # target: symptom time
                becomes_infected=rep(-Inf,n), # when node becomes infectious
                infected_by=rep(-Inf,n), # record for later
                done_infecting=rep(F,n)) # flag for propagation for current node is done

# init source
nodes$infectious[[s]] <- T
nodes$incubating_until[[s]] <- rexp(1,1/beta_i[[s]])

# start sim
resolution <- 1 # how smoothly will the propagation occur
for(t in seq(1,t_limit,resolution)){
  # get list of currently infectious
  infected <- nodes %>% filter(infectious & !done_infecting)
  
  # loop through them and set infectious times for all neighbors
  for(node in infected$node){
    # set flag
    nodes$done_infecting[[node]] <- T
    
    # get neighbors
    neighbors <- which(G_M[node,]==1)
    
    # loop through and infect neighbors
    for(neighbor in neighbors){
      if(!nodes$infectious[[neighbor]]){ # no loops!
        nodes$becomes_infected[[neighbor]] <- rexp(1,1/alpha_ij[node,neighbor]) + t
        nodes$infected_by[[neighbor]] <- node
      }
    }
  }
  
  # check to see if any nodes become infected
  check <- nodes %>% filter(!infectious & (becomes_infected > 0))
  for(node in check$node){
    # check times
    if(abs(nodes$becomes_infected[[node]] - t) <= resolution){
      # if true, infect and set incubation time
      nodes$infectious[[node]] <- T
      nodes$incubating_until[[node]] <- rexp(1,1/beta_i[[node]]) + t
    }
  }
}

# pull valid infection times (if not infected, set time to Inf)
infection_times <- nodes %>% select(node,incubating_until) %>%
  mutate(incubating_until=ifelse((0 < incubating_until) & (incubating_until < t_limit),
                                 incubating_until,Inf) ) %>%
  rename(infection_time = incubating_until)

# save for later
true_s <- s
save(infection_times,G,beta_i,alpha_ij,true_s,file="sias_spread42")
