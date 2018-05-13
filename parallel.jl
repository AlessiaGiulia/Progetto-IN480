@everywhere function pt(args...)
	d=length(args)
	mat=SharedArray(eye(d+1))
	@parallel for i in range(1,d)
		mat[i,d+1]=args[i]
	end
	return mat
end

#____________________________________________________________________________________________________________________________

@everywhere function ps(args...)
	d=length(args)
	mat=SharedArray(eye(d+1))
	@parallel for k in range(1,d)
		mat[k,k]=args[k]
	end
	return mat
end

#____________________________________________________________________________________________________________________________


@everywhere function pr(args...)
    args = collect(args)
    n = length(args)
	mat=eye(3) #vedere se farla diventare una sharedArray
    if n == 1 # rotation in 2D
	@async begin
		angle = args[1]
		#COS=remotecall_fetch(cos,2,angle)
			#COS=remotecall(cos,2,angle)
		COS = cos(angle)
		SIN = sin(angle)
		@sync begin
			
        		mat[1,1] = COS
			mat[1,2] = -SIN
        		mat[2,1] = SIN
			mat[2,2] = COS
		end
	
   end
   end
	if n == 3 # rotation in 3D
        	@async begin
			angle = norm(args); axis = normalize(args)
        		COS = cos(angle); SIN= sin(angle)
			@sync begin
        		if axis[2]==axis[3]==0.0    # rotation about x
         		   mat[2,2] = COS;    mat[2,3] = -SIN;
          		   mat[3,2] = SIN;    mat[3,3] = COS;

       			elseif axis[1]==axis[3]==0.0   # rotation about y
           			mat[1,1] = COS;    mat[1,3] = SIN;
           			mat[3,1] = -SIN;    mat[3,3] = COS;
       			elseif axis[1]==axis[2]==0.0    # rotation about z
            			mat[1,1] = SIN;    mat[1,2] = -SIN;
            			mat[2,1] = COS;    mat[2,2] = COS;
        
        		else
	   			I=eye(3); u=axis
	   			Ux=[0 -u[3] u[2] ; u[3] 0 -u[1] ;  -u[2] u[1] 1]
	   			UU =[u[1]*u[1]    u[1]*u[2]   u[1]*u[3];
             			     u[2]*u[1]    u[2]*u[2]   u[2]*u[3];
            			     u[3]*u[1]    u[3]*u[2]   u[3]*u[3]]
	    			mat[1:3,1:3]=COS*I+SIN*Ux+(1.0-COS)*UU
			end
			end
		end
	end
	return mat
end

#____________________________________________________________________________________________________________________________

function premoveDups(CW)
	CW=collect(Set(CW))
	CWs=collect(pmap(sort,CW))
	no_duplicates=Dict()
	@parallel for f in CWs
		no_duplicates[f] = []
	end
	@parallel for f in CW
		no_duplicates[sort(f)]=[f]
	end

	@parallel for f in values(no_duplicates)
	append!(CW,f[1])
	end
	return CW
end

#____________________________________________________________________________________________________________________________

 
@everywhere function larRemoveVertices(V,FV)
	vertDict= Dict()
	index,defaultValue,CW,W,FW = -1,-1,[],[],[]
	for (k,incell) in enumerate(FV)
		outcell=[]
		for v in incell
			key=vcode(4)(V[v+1])
			if get(vertDict,key,defaultValue)==defaultValue
				index =index+1
				vertDict[key]=index
				outcell =append!(outcell,index)
				W= append!(W,[eval(parse(key))])                   
			else
				outcell= append!(outcell,vertDict[key])
			end
		end
		FW =append!(FW,[outcell])
	end
	return W,FW
end