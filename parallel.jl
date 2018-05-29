function pprepKey(args)
               v=join(args,",")
return(v)
end


#____________________________________________________________________________________________________________________________

function pfixedPrec(PRECISION)
	function pfixedPrec0(value) 
		out=round.(value,PRECISION)
		if out==-0.0
			out=0.0
		end
		return string(out)
	end
	return pfixedPrec0
end

#____________________________________________________________________________________________________________________________

function pvcode(PRECISION=4)
	function pvcode0(vect)
		#return prepKey(pmap(fixedPrec(PRECISION),vect)) dovrebbe essere così ma in julia non serve mettere prepKey per ottenere lo stesso risultato di python probabilmente perchè il map funziona in maniera differente
		return fixedPrec(PRECISION)(vect) #quando vect è contiene più array bisogna usare map
	end
	return pvcode0
end

#____________________________________________________________________________________________________________________________


function t(args...)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
        	mat[k,d+1]=args[k]
	end
	return mat
end

#____________________________________________________________________________________________________________________________

function s(args...)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
		mat[k,k]=args[k]
	end
	return mat
end


#____________________________________________________________________________________________________________________________


function r(args...)
    args = collect(args)
    n = length(args)

    if n == 1 # rotation in 2D
        angle = args[1]; COS = cos(angle); SIN = sin(angle)
        mat = eye(3)
        mat[1,1] = COS;    mat[1,2] = -SIN;
        mat[2,1] = SIN;    mat[2,2] = COS;
    end

     if n == 3 # rotation in 3D
        mat = eye(4)
        angle = norm(args); axis = normalize(args)
        COS = cos(angle); SIN= sin(angle)
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

return mat

end


#____________________________________________________________________________________________________________________________

function premoveDups(CW)
	CW=collect(Set(CW))

	CWs=collect(@sync pmap(sort,CW))
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

 
@everywhere function plarRemoveVertices(V,FV)
	vertDict= Dict()
	index,defaultValue,CW,W,FW = -1,-1,[],[],[]
	@async begin
	       for (k,incell) in enumerate(FV)
	       	   	outcell=[]
			@sync begin
			      for v in incell
			      	  key=pvcode(4)(V[v+1])
					if get(vertDict,key,defaultValue)==defaultValue
				   	   index =index+1
				   	   vertDict[key]=index
				   	   append!(outcell,index)
				   	   append!(W,[eval(parse(key))])
					else
						append!(outcell,vertDict[key])
					end
				end
			end
			append!(FW,[outcell])
		end
	end
	return W,FW
end


#____________________________________________________________________________________________________________________________


#classe pStruct



@everywhere type pStruct
	body::Array
	box
	name::AbstractString
	dim
	category::AbstractString
	
	function pStruct()
		self=new([],Nullable{Any},"new",Nullable{Any},"feature")
		self.name=string(object_id(self))
		return self

	end

	function pStruct(data::Array)
		self=pStruct()
		self.body=data
		#println(pbox(data))
		self.box=pbox(self)
		self.dim=length(self.box[1])
		return self
	end
	
	function pStruct(data::Array,name)
		self=pStruct()
		self.body=[item for item in data]
		self.box=pbox(self)
		self.dim=length(self.box[1])
		self.name=string(name)
		return self
	end

	function pStruct(data::Array,name,category)
		self=pStruct()
		self.body=[item for item in data]
		self.box=pbox(self)
		self.dim=length(self.box[1])
		self.name=string(name)
		self.category=string(category)
		return self
	end
	
end

	function name(self::pStruct)
		return self.name
	end
	function category(self::pStruct)
		return self.category
	end
	
	function len(self::pStruct)
		return length(self.body)
	end
	function getitem(self::pStruct,i::Int)
		return self.body[i]
	end
	function setitem(self::pStruct,i,value)
		self.body[i]=value
	end
	function pprint(self::pStruct)
		return "<Struct name: $(self.__name__())"
	end
	function set_name(self::pStruct,name)
		self.name=string(name)
	end
	function clone(self::pStruct,i=0)
		newObj=deepcopy(self)
		if i!=0
			newObj.name="$(self.__name__())_$(string(i))"
		end
		return newObj
	end
	function set_category(self::pStruct,category)
		self.category=string(category)
	end

#____________________________________________________________________________________________________________________________

@everywhere function plarBoundary(self)
	data=pstruct2lar(self)
	if length(data)==3
		V,FV,EV=data
		return V,FV,EV
	else
		return "<Struct name $(self.__name__()) : boundary non computable"
	end
end


#____________________________________________________________________________________________________________________________


@everywhere function pstruct2lar(structure)
	listOfModels=pevalStruct(structure)
	vertDict= Dict()
	index,defaultValue,CW,W,FW = -1,-1,[],[],[]
	
	for model in listOfModels
		if  length(model)==2
			V,FV=model
		elseif lenght(model)==3
			V,FV,EV=model
		end
		@sync begin
			for (k,incell) in enumerate(FV)
				outcell=[]
			@async begin
			      for v in incell
			      	  key=pvcode(4)(V[v+1])
				  if get(vertDict,key,defaultValue)==defaultValue
					index =index+1
                    			vertDict[key]=index
					append!(outcell,index)
					append!(W,[eval(parse(key))])
				  else
					append!(outcell,vertDict[key])
				  end
			       end	
			end
			append!(CW,[outcell])
			end
		end
		if length(model)==3
		   	@sync begin
				for (k,incell) in enumerate(FV)
					outcell=[]
					@async begin
				       	       for v in incell
				       	       	   key=pvcode(4)(V[v+1])
							if get(vertDict,key,defaultValue)==defaultValue
							   index =index+1
							   vertDict[key]=index
							   append!(outcell,[index])
							   append!(W,[eval(parse(key))])                   
					   		else
								append!(outcell,vertDict[key])
							end
						end
					end
				append!(FW,[outcell])
				end
			end
		end
	end
	
	if length(listOfModels[end])==2
		if length(CW[1])==2
			CW=pmap(Tuple,pmap(sort,CW))
		else
			CW=premoveDups(CW)
		end
		return W,CW
	end
	
	if length(listOfModels[end])==3
		FW=pmap(Tuple,pmap(sort,FW)) #da controllare nel test
		CW=premoveDups(CW)
		return W,CW,FW
	end
end



#____________________________________________________________________________________________________________________________


@everywhere function pembedTraversal(cloned,obj,n,suffix)

	for i in range(1,len(obj))
	      	 if (isa(obj.body[i],Matrix) || isa(obj.body[i],SharedArray))
			mat=obj.body[i]
			d,d=size(mat)
			newMat=eye(d+n*1)
			@sync begin
			for h in range(1,d-1)
				@async begin
				       for k in range(1,d-1)
					   newMat[h,k]=mat[h,k]
					end
				end
				newMat[h,d-1+n*1]=mat[h,d-1]
			end
			end
			append!(cloned.body,newMat)


		elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&& length(obj.body[i])==2 
			V,EV=deepcopy(obj.body[i])
			dimadd=fill([0.0],n)
			@sync begin
			      for k in dimadd
			      	  @async begin
				       for v in V
				       	   append!(v,k)
					end
				end
				end
			end
			append!(cloned.body,[(V,EV)])
		elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&& length(obj.body[i])==3 #potrebbe dare problemi in ndarray
			V,FV,EV=obj.body[i]
			dimadd=fill([0.0],n)
			@sync begin
			      for k in dimadd
			      	  @async begin
				       for v in V
				       	   append!(v,k)
					end
				end
				end
			end
			append!(cloned.body,[(V,FV,EV)])
		 
		elseif isa(obj.body[i],pStruct)
			newObj=pStruct()
			@async begin
			newObj.box=hcat((obj.body[i].box,[fill([0],n),fill([0],n)]))
			newObj.category=obj.body[i].category
			append!(cloned.body,embedTraversal(newObj,obj.body[i],n,suffix))
			end
		end
	end
	return cloned
end

#____________________________________________________________________________________________________________________________


@everywhere function pembedStruct(n)
	function pembedStruct0(self,suffix="New")
		if n==0
			return self, length(self.box[1])
		end
		cloned=pStruct()
		cloned.box=hcat((self.box,[fill([0],n),fill([0],n)]))	
		cloned.name=string(object_id(cloned))
		cloned.category=self.category
		cloned.dim=self.dim+n
		cloned=pembedTraversal(cloned,self,n,suffix)
		return cloned
	end
	return pembedStruct0
end

#____________________________________________________________________________________________________________________________

function pbox(model)
	if isa(model,Matrix)
		return []
	elseif isa(model,pStruct)
		dummyModel=deepcopy(model)
		dummyModel.body=Any[]
		for term in model.body 
			if isa(term,pStruct)
				push!(dummyModel.body,[term.box,[0,1]])
			else
				push!(dummyModel.body,term)
			end
		end
		listOfModels=pevalStruct(dummyModel)
		#dim=checkStruct(listOfModels)
		theMin,theMax=pbox(listOfModels[1])
		
		for theModel in listOfModels[2:end]
			modelMin,modelMax= pbox(theModel)
			
			for (k,val) in enumerate(modelMin)
				#
				if (val < theMin[k])
					theMin[k]=val
				end
			end
			for (k,val) in enumerate(modelMax)
				if (val > theMax[k])
					theMax[k]=val
				end
			end
		end
		return Array[theMin,theMax]

	elseif (isa(model,Tuple) ||isa(model,Array))&& (length(model)==2 || length(model)==3)
		V=model[1]
	theMin=[]
	theMax=[]
		for j in range(1,length(V[1]))
			Min=V[1][j]
			Max=V[1][j]
			for i in range(1,length(V))
				Min=min(Min,V[i][j])
				Max=max(Max,V[i][j])
			end
			push!(theMin,Min)
			push!(theMax,Max)
		end

	return Array[theMin,theMax]
	end
end

 

#____________________________________________________________________________________________________________________________


@everywhere function plarApply(affineMatrix)
	function plarApply0(model)
	    if length(model)==2
	       V,CV=deepcopy(model)
	    elseif length(model)==3
      	    	   V,CV,FV = deepcopy(model)
	    end
	   V1=Array{Float64}[]
	   V1= @sync @parallel (append!)for v in V
	    	       append!(v,[1.0])
	    	      [collect(vec((v')*transpose(affineMatrix)))]
		end
		
		 
		 pmap(x->pop!(x),V1)
		   
		  
	     
	     if length(model)==2
		return V1,CV
	     elseif length(model)==3
		return V1,CV,FV
	      end
	 end 
	 return plarApply0
end

#____________________________________________________________________________________________________________________________

@everywhere function pcheckStruct(lst)
	obj = lst[1]
	if (isa(obj,Matrix) || isa(obj,SharedArray))
		dim=size(obj)[1]-1
	elseif(isa(obj,Tuple) || isa(obj,Array))
		dim=length(obj[1][1])
	
	elseif isa(obj,pStruct)
		dim=length(obj.box[1])
	end
	return dim
end		



#____________________________________________________________________________________________________________________________

@everywhere function ptraversal(CTM,stack,obj,scene=[])
	for i in range(1,len(obj))
		if (isa(obj.body[i],Matrix)|| isa(obj,SharedArray))
			CTM=CTM*obj.body[i]
		elseif (isa(obj.body[i],Tuple) || isa(obj.body[i],Array)) && (length(obj.body[i])==2 || length(obj.body[i])==3)
			l=plarApply(CTM)(obj.body[i])
			push!(scene,l)
		elseif isa(obj.body[i],pStruct)
			push!(stack,CTM)	
			traversal(CTM,stack,obj.body[i],scene)
			CTM=pop!(stack)
			
		end
	end
	return scene
end


#____________________________________________________________________________________________________________________________
@everywhere function pevalStruct(self)
	dim = pcheckStruct(self.body)
   	CTM, stack = eye(dim + 1), []
   	scene = ptraversal(CTM, stack, self, []) 
return scene
end


