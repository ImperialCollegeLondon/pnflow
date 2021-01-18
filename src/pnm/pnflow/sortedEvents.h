#ifndef SORTEDEVENTS_H
#define SORTEDEVENTS_H


/**
// The class takes care of all the sorted events
*/
template<typename Type>
class SortedEventsBase
{
public:

	void quickInsert(Type elm) {sortedContainer_.push_back(elm);}
	virtual void insert(Type elm) = 0;
	virtual bool remove(Type elm) = 0;
	virtual void sortEvents() = 0;
	Type pop()	{
		ensure(!sortedContainer_.empty());
		Type elem = sortedContainer_.back();
		sortedContainer_.pop_back();
		return elem;
	};
	Type peek() const	{
		ensure(!sortedContainer_.empty());
		return sortedContainer_.back();
	};
	bool present(Type elem) const	{
		for(size_t i = 0; i < sortedContainer_.size(); ++i)
			if(elem == sortedContainer_[i])
				return true;
		return false;
	};
	Type at(size_t pos) const	{
		ensure(sortedContainer_.size() > pos);
		return sortedContainer_[pos];
	};
	size_t size() const {return sortedContainer_.size();}
	bool empty() const {return sortedContainer_.empty();}
	void clear() {sortedContainer_.clear();}

	bool checkIfThere(Type elem) const	{    // Very expensive function that might however be useful
		for(size_t i = 0; i < sortedContainer_.size(); ++i)
			if(elem == sortedContainer_[i]) return true;
		return false;
	};
protected:

													                  // during debugging


	typedef typename std::vector< Type >::iterator        EventItr;

	std::vector<Type>                            sortedContainer_;
};




/**
// The class takes care of all the sorted events
*/
template<typename Type, typename CFunc>
class SortedEvents  : public SortedEventsBase<Type>
{
public:

	void insert(Type elem){
		//ensure(!checkIfThere(elem));    // Only enable for emergency debugging  => EXTREMELY expensive
		 //if (checkIfThere(elem))
		 //{ cout<<"reins"<<elem->isInWatFloodVec()<<" ";cout.flush(); //return;
			//// ((Type*)(&elem))[1000000000]=0;
			//remove(elem);
		 //}
		SortedEventsBase<Type>::sortedContainer_.insert(
			lower_bound(SortedEventsBase<Type>::sortedContainer_.begin(),
			SortedEventsBase<Type>::sortedContainer_.end(),
			elem, CFunc()), elem);
	};
	bool remove(Type elem)
	{
		EventItr delCand = lower_bound(SortedEventsBase<Type>::sortedContainer_.begin(), SortedEventsBase<Type>::sortedContainer_.end(), elem, CFunc());
		while(delCand != SortedEventsBase<Type>::sortedContainer_.end() && *delCand != elem) ++delCand;
		if (delCand != SortedEventsBase<Type>::sortedContainer_.end())
		{
			SortedEventsBase<Type>::sortedContainer_.erase(delCand);
			return true;//SortedEventsBase<Type>::sortedContainer_.empty();
		}
		//else if (delCand != SortedEventsBase<Type>::sortedContainer_.end() && *delCand != elem) 
		//{
			//std::cout<<"  Error:*delCand==elem,"<<checkIfThere(elem)<<"  ";std::cout.flush();
			//return false;
		//}

		return false;
	};
	void sortEvents(){
	sort(SortedEventsBase<Type>::sortedContainer_.begin(), SortedEventsBase<Type>::sortedContainer_.end(), CFunc());
};


private:

	typedef typename std::vector< Type >::iterator        EventItr;

};



#endif
