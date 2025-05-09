module scheduler 

---------- Filesystem ----------------
    sig File {}

    one sig Filesystem {
        files: seq File,
    }
    {
        File in files.elems
    }
---------------------------------------

abstract sig State {}
one sig NE extends State {}
one sig S extends State {}
one sig C extends State {
    var commit_order : seq Command
}
one sig W extends State {}
one sig E extends State{}
one sig D extends State{} 

// Just check that commit order respects actual partial Order


sig Command {
   
     var command_state: one State,
   
    -- This is the syntactic order predicted by the preprocessor. 
    syntactic_order : set Command,    


    // Details the operation_on_filesystems of each command from File to File for each possible state.
    // A read is of the form of an identity transformation for a File (ie content is unchanged)
    // while a write involves changing the content of a file!

    // Files Command reads from 
    var read_set : set File ,

    // Files Command writes from 
    var write_set : set File  ,

    // var commit_order : lone Command ,

    var side_effect : lone Int , 

    var execute_idx : lone Int

}
----------------------------------------------------------------------------------------------------
// Well Formedness

    pred partialOrder[r: set (Command -> Command)] {
        no r & iden // Anti-reflexivity
        all x: Command, y: Command, z: Command |
        {
            (x->y in r and y->z in r) implies x->z in r // Transitivity
        }
        no (r & ~r) // anti-symmetric
    }

    pred wellFormed {
        partialOrder[syntactic_order]

        not (C.commit_order.isEmpty) implies C.commit_order'.subseq[0,C.commit_order.lastIdx] = C.commit_order
        
        all c : Command | {
            c in (C.commit_order'.elems - C.commit_order.elems) iff {
                (commit_node[c] or direct_commit[c])
            }
        }
        #C.commit_order'.elems = add[#C.commit_order.elems ,#(C.commit_order'.elems - C.commit_order.elems)]                     
    }

    
------------------------------------------------------------------------------------------------------

// // Helpers
//     cmd addActionOrder(c : Command) { 

//     }
    fun command_pred[c : Command] : set Command{
        {x : Command | 
            c in x.^syntactic_order
        }
    }
    fun Frontier_set [ordr : Command -> Command] : set Command {

        {x : Command | 
            {
                all y : Command | x in y.^ordr implies {committed[y]}
            }
        }
    }



    pred hasDependency [c1 : Command , c2 : Command ] { 
        // Forward Dependency C1 writes , C2 reads
        some (c1.write_set & c2.read_set) or 
        // Backward Dependency C1 reads,  C2 writes
        some (c1.read_set & c2.write_set) or 
        // Write dependency C1 writes, C2 writes
        some (c1.write_set & c2.write_set)
    }

    pred maintainRW[c : Command] {
        c.read_set' = c.read_set 
        c.write_set' = c.write_set
    }

    fun nonCommittedDependencies[c : Command] : set Command
    {
        {x : Command | (c in x.^syntactic_order) and !committed[x] and hasDependency[x,c] }
    }

    pred setExcIdx [c : Command ] {
        C.commit_order.isEmpty implies c.execute_idx' = 0 
        not C.commit_order.isEmpty implies c.execute_idx' = C.commit_order.lastIdx 
    }
    pred action_order_violated[c : Command] {
           
        some c2 : Command  | { 
                c2 in C.commit_order.elems
                c2 in command_pred[c]
                hasDependency[c2,c]
                c.execute_idx < C.commit_order.idxOf[c2]
            }

    }

   
---------------------------------------------------------------------------------------------------------
//State transition rules


    // NE -> NE 
    pred stayNE[c : Command] { 
        c.command_state = NE
        after (c.command_state = NE)
        
        c not in Frontier_set[syntactic_order]
       
    }
    // NE -> D
    pred direct_execute [c : Command] {
        c.command_state = NE
        after (c.command_state = D)
        c in Frontier_set[syntactic_order]
        some c.side_effect
        setExcIdx[c]
    }

    // D -> D
    pred direct_executing [c : Command] {
        c.command_state = D
        after (c.command_state = D)
        eventually  (c.command_state != D)
    }

    pred direct_commit [c : Command] {
        c.command_state = D
        after(c.command_state = C)
        
    }

    // NE -> E 
    pred execute_command[c : Command] {
        c.command_state = NE
        after (c.command_state = E)
        no c.side_effect
        // TODO : what if a command gets commited at the same time . This is probably correct
        // as we are being extra cautious
         setExcIdx[c]
        
    }

    // E -> NE 
    pred side_effect_found[c : Command] {
        c.command_state = E
        after (c.command_state = NE)
        some c.side_effect
     }

    // E -> E
    pred stayE [c : Command] { 
        c.command_state = E 
        after(c.command_state = E)
        eventually(c.command_state != E)
        no c.side_effect
    }
    // E -> W
    pred command_exec_finished [c : Command] { 
        c.command_state = E
        after(c.command_state = W)
        no c.side_effect
    }

    // W -> W 
    pred command_waiting [c : Command] {
        c.command_state = W
        after(c.command_state = W)
        
        some c2 : Command | {
            c in c2.^syntactic_order
            (c2.command_state = NE or c2.command_state = E)
        }

    }
    // W -> NE 
    pred trace_executor_found_dependency[c : Command] {
        c.command_state = W
        after (c.command_state = NE)
        all c2 : Command | {
            c in c2.^syntactic_order
            (c2.command_state != NE or c2.command_state != E)
        }
        some nonCommittedDependencies[c]
    }
    // W -> S
    pred mark_speculated[c : Command] {
        c.command_state = W
        after (c.command_state = S)
        
        no nonCommittedDependencies[c]

    }

    // S -> C
    pred commit_node[c : Command] {
        c.command_state = S
        after (c.command_state = C )

        c in Frontier_set[syntactic_order]
        not action_order_violated[c]

    }

    // S -> S
    pred staySpec[c : Command] {
        c.command_state = S
        after (c.command_state = S )

        c not in Frontier_set[syntactic_order]

    }
    // S -> NE
    pred speculated_dep_found[c : Command] {
        c.command_state = S
        after (c.command_state = NE )
        action_order_violated[c]
    }
    pred NEaction[c : Command]  {
        stayNE[c] or execute_command[c] or direct_execute[c]
    }

    pred Daction [c : Command] {
        direct_executing[c] or direct_commit[c]
    }

    pred Eaction [c : Command] { 
        stayE[c] or command_exec_finished[c] or side_effect_found[c]
    }
    pred Waction [c : Command] {
        command_waiting[c] 
        or trace_executor_found_dependency[c] 
        or mark_speculated[c]
    } 
    pred Saction [c : Command] {
        commit_node[c] 
        or speculated_dep_found[c] 
        or staySpec[c]
    }

    pred committed[c: Command] {
        (c.command_state) = C 
        {(c.command_state) = C} implies always (c.command_state = C)  
    }

    pred validAction[c : Command] {
       NEaction[c] or Daction[c]  or Eaction[c] or Waction[c]  or Saction[c]   or committed[c]
       not maintainRW[c] implies (execute_command[c] or direct_execute[c])
       not(c.side_effect' = c.side_effect) implies execute_command[c]
       not (c.execute_idx' = c.execute_idx) implies (execute_command[c] or direct_execute[c])
       
    }
------------------------------------------------------------------------------------------------

// Scheduler Behavior
    pred traces {
        always wellFormed
        all c : Command | {
             always validAction[c]
        }
    }

    // Initial state
    pred init {
        all c : Command | {
            c.command_state = NE
            no c.read_set 
            no c.write_set
            no c.side_effect
            no c.execute_idx
        }
        
        #C.commit_order.elems = 0
    }

    // Scheduler behavior
    pred scheduler_e2e {
        init
        traces
    }

    // All commands have been committed at the end
    pred final {
        all c : Command | committed[c]
        
    }