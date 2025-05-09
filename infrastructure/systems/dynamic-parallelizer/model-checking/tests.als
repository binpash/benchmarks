
open scheduler as sch
open util/ordering[Command] as lin
run satrun { 
    scheduler_e2e
} for exactly 6 State, exactly 6 Command, 6 File,6 seq

run sideffect_allowed {
    scheduler_e2e 
    eventually(some c : Command | some c.side_effect)
    }  for exactly 6 State, 2 Command, 6 File ,2 seq


run no_syntactic_order { 
    scheduler_e2e
    no syntactic_order
    
} for exactly 6 State, exactly 6 Command, 6 File ,6 seq

run allow_io {
    scheduler_e2e
    eventually (some read_set)
    eventually (some write_set)
} for exactly 6 State, exactly 6 Command, 6 File,6 seq

run change_write_set {
    scheduler_e2e
    some f : File , c : Command {
        eventually(f in c.write_set until f not in c.write_set)
    }
} for exactly 6 State, exactly 6 Command, 6 File,6 seq

   pred commit_precede[c1 : Command,c2 : Command] {
        ((c1 + c2) in C.commit_order.elems and 
        C.commit_order.idxOf[c1] < C.commit_order.idxOf[c2])
    }
pred extra_preservation{ 
        all disj c1,c2 : Command | {
            commit_precede[c1,c2] implies not hasDependency[c2,c1] 
        }
    }
run no_extra_preservation { scheduler_e2e 
     eventually(not extra_preservation)} for exactly 6 State, 6 Command, 6 File,6 seq
