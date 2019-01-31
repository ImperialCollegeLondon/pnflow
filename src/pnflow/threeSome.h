#ifndef THREESOME_H
#define THREESOME_H


/**
// A little storage class for three elements
*/
template<typename TOne, typename TTwo, typename TThree>
class ThreeSome
{
public:

    ThreeSome() {}
    ThreeSome(TOne one, TTwo two, TThree three) : m_first(one), m_second(two), m_third(three) {}

    const TOne& first() const {return m_first;}
    const TTwo& second() const {return m_second;}
    const TThree& third() const {return m_third;}
    void first(TOne entry) {m_first = entry;}
    void second(TTwo entry) {m_second = entry;}
    void third(TThree entry) {m_third = entry;}
	bool operator==(ThreeSome b)
	{return m_first == b.first() && m_second == b.second() && m_third == b.third() ;}
	bool operator!=(ThreeSome b)
	{return m_first != b.first() || m_second != b.second() || m_third != b.third() ;}

private:

    TOne                m_first;
    TTwo                m_second;
    TThree              m_third;

};


/**
// A little storage class for four eleemnts
*/
template<typename TOne, typename TTwo, typename TThree, typename TFour>
class FourSome
{
public:

    FourSome() {}
    FourSome(TOne one, TTwo two, TThree three, TFour four) : m_first(one), m_second(two), m_third(three), m_fourth(four) {}

    const TOne& first() const {return m_first;}
    const TTwo& second() const {return m_second;}
    const TThree& third() const {return m_third;}
    const TFour& fourth() const {return m_fourth;}
    void first(TOne entry) {m_first = entry;}
    void second(TTwo entry) {m_second = entry;}
    void third(TThree entry) {m_third = entry;}
    void fourth(TFour entry) {m_fourth = entry;}

private:

    TOne                m_first;
    TTwo                m_second;
    TThree              m_third;
    TFour               m_fourth;

};

#endif
