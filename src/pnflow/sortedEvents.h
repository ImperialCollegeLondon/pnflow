#ifndef SORTEDEVENTS_H
#define SORTEDEVENTS_H


/**
// The class takes care of all the sorted events
*/
template<typename Type>
class SortedEventsBase
{
public:

    void quickInsert(Type elm) {m_sortedContainer.push_back(elm);}
    virtual void insert(Type elm) = 0;
    virtual bool remove(Type elm) = 0;
    virtual void sortEvents() = 0;
    Type pop()	{
		softAssert(!m_sortedContainer.empty());
		Type elem = m_sortedContainer.back();
		m_sortedContainer.pop_back();
		return elem;
	};
    Type peek() const	{
		softAssert(!m_sortedContainer.empty());
		return m_sortedContainer.back();
	};
    bool present(Type elem) const	{
		for(size_t i = 0; i < m_sortedContainer.size(); ++i)
			if(elem == m_sortedContainer[i])
				return true;
		return false;
	};
    Type at(size_t pos) const	{
		softAssert(m_sortedContainer.size() > pos);
		return m_sortedContainer[pos];
	};
    size_t size() const {return m_sortedContainer.size();}
    bool empty() const {return m_sortedContainer.empty();}
    void clear() {m_sortedContainer.clear();}

    bool checkIfThere(Type elem) const	{    // Very expensive function that might however be useful
		for(size_t i = 0; i < m_sortedContainer.size(); ++i)
			if(elem == m_sortedContainer[i]) return true;
		return false;
	};
protected:

                                                                      // during debugging

	
    typedef typename std::vector< Type >::iterator        EventItr;

    std::vector<Type>                            m_sortedContainer;
};




/**
// The class takes care of all the sorted events
*/
template<typename Type, typename CFunc>
class SortedEvents  : public SortedEventsBase<Type>
{
public:

    void insert(Type elem){
		//softAssert(!checkIfThere(elem));    // Only enable for emergency debugging  => EXTREMELY expensive
		 //if (checkIfThere(elem))
		 //{ cout<<"reins"<<elem->isInWatFloodVec()<<" ";cout.flush(); //return;
			//// ((Type*)(&elem))[1000000000]=0;
			//remove(elem);
		 //}
		SortedEventsBase<Type>::m_sortedContainer.insert(
			lower_bound(SortedEventsBase<Type>::m_sortedContainer.begin(),
			SortedEventsBase<Type>::m_sortedContainer.end(),
			elem, CFunc()), elem);
	};
    bool remove(Type elem)
    {
		EventItr delCand = lower_bound(SortedEventsBase<Type>::m_sortedContainer.begin(), SortedEventsBase<Type>::m_sortedContainer.end(), elem, CFunc());
		while(delCand != SortedEventsBase<Type>::m_sortedContainer.end() && *delCand != elem) ++delCand;
		if (delCand != SortedEventsBase<Type>::m_sortedContainer.end())
		{
			SortedEventsBase<Type>::m_sortedContainer.erase(delCand);
			return true;//SortedEventsBase<Type>::m_sortedContainer.empty();
		}
		//else if (delCand != SortedEventsBase<Type>::m_sortedContainer.end() && *delCand != elem) 
		//{
			//std::cout<<"  Error:*delCand==elem,"<<checkIfThere(elem)<<"  ";std::cout.flush();
			//return false;
		//}

		return false;
	};
    void sortEvents(){
    sort(SortedEventsBase<Type>::m_sortedContainer.begin(), SortedEventsBase<Type>::m_sortedContainer.end(), CFunc());
};


private:

    typedef typename std::vector< Type >::iterator        EventItr;

};



#endif
