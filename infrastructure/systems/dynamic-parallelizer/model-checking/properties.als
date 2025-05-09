
open scheduler as sch
open util/ordering[Command] as lin

// Safety:


    //If a command has been committed, then all of its dependencies have also been committed.
    pred commit_precede[c1 : Command,c2 : Command] {
        ((c1 + c2) in C.commit_order.elems and 
        C.commit_order.idxOf[c1] < C.commit_order.idxOf[c2])
    }
    
    pred dependency_preservation {
        all c1,c2 : Command | {
            commit_precede[c1,c2] implies not(c1 in c2.^syntactic_order)
        }
    }
    pred commit_order_contained {
        some c : Command { 
            (c in C.commit_order.elems ) iff committed[c]
        }
        final implies C.commit_order.elems = Command
    }
    pred action_order_maintained { 
        no c1,c2 : Command  | { 
                (c1 + c2) in C. commit_order.elems
                c2 in c1.^syntactic_order // c1->c2
                hasDependency[c1,c2]
                c2.execute_idx < C.commit_order.idxOf[c1]
            }
        
    }

    check action_order_maintained { scheduler_e2e implies always(action_order_maintained)} for exactly 6 State, 6 Command, 6 File,6 seq
    check dep_preserve { scheduler_e2e implies always(dependency_preservation)} for exactly 6 State, 6 Command, 6 File,6 seq 
    
// Termination

    // Once terminated, nothing is scheduled.
    check termination_is_final {final implies (always final)  } for exactly 6 State, 6 Command, 6 File , 6 seq 
    // This shows that the scheduler terminates, with all Commands committed.
    check sched_terminates {scheduler_e2e implies (eventually final) } for exactly 6 State, 6 Command, 6 File, 6 seq 

    check commit_order_contained {scheduler_e2e implies always commit_order_contained} for exactly 6 State, 6 Command, 6 File, 6 seq 

    // Only 1 command should be executing outside the sandbox at a time
    // Currently failing - there are fixes but not implemented as afaik the code doesn't handle this yet
    // check {scheduler_e2e implies always(lone c : Command | Daction[c])}  for exactly 6 State, 6 Command, 6 File, 6 seq 

    // TODO : Dependency between commands on the frontier

