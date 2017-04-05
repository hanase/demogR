\name{leslie.matrix}
\alias{leslie.matrix}

\title{ Leslie Matrix }
\description{
  Generates a Leslie matrix for demographic projection from vectors of
  age-specific cumulative survival and fertility.
}
\usage{
leslie.matrix(lx, mx, 
		    L = TRUE, 
		    peryear = 5, 
		    one.sex = TRUE, 
		    SRB = 1.05, 
		    infant.class = TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lx}{ vector of either age-specific cumulative survival or
  person-years lived in the interval }
  \item{mx}{ age-specific fertility rates }
  \item{L}{ logical; if 'TRUE', \code{lx} is taken to be person-years
  lived in the interval nLx, while if 'FALSE', \code{lx} is taken to
  be cumulative survival to exact age x+n.  Defaults to 'TRUE'. }
  \item{peryear}{ Multiplier for fertility.  Defaults to \code{peryear=5}. }
  \item{one.sex}{ logical; If 'TRUE', fertility rates will be divided
  by 1/(1+SRB). }
  \item{SRB}{ sex ratio at birth.  Defaults to \code{SRB=1.05}. }
  \item{infant.class}{ logical;  'TRUE' if \code{lx} contains a value
  for the infant age-class.   }
}

\details{ 

  Constructs a k x k age-classified demographic projection matrix with
  age-specific survival probabilities along the sub-diagonal and
  age-specific fertilities along the first row of the
  matrix. 

  \code{lx} and \code{mx} are assumed to be of the same
  length.  The resulting matrix is truncated to insure that there are
  no post-reproductive classes.  This is important for ensuring
  irreducibility of the resulting matrix.

  If \code{mx} is longer than \code{lx}, \code{mx} is trucated to be
  the same length as \code{lx}.  If \code{lx} is longer than
  \code{mx}, a warning is issed and \code{lx} is truncated to be the
  same length as \code{mx}.

  Fertility is assumed to be birth-flow (Caswell 2001).  That is,
  breeding is assumed to be continuous and the individual elements of
  the first row of the Leslie matrix are averaged over successive
  age-classes.  Fertility rates are typically given in annualized
  form.  If this is the case and the age-classes are wider than one
  year, then \code{peryear} can be used to appropriately scale up the
  annual values.

  The default behavior is to use person-years lived in the interval as
  the survival measure.   If \code{infant.class=TRUE}, lx is taken to
  have a value for the infant age class (i.e., a shorter class width
  than the other elements of \code{lx}.  What is done when there is an
  infant class depends on what the values in \code{lx} represent. If
  \code{L=TRUE}, then the first two values of lx are combined to form
  the total person-years for the first age-class in the Leslie
  matrix.  Human demographic data from abridged life tables typically
  come with age classes x = 0, 1, 5, 10, ... Thus, combining the
  person-years for the first two age classes gives an initial age
  class of the correct width. If \code{infant.class=TRUE} and
  \code{L=FALSE}, the second element of \code{lx} is deleted.
  Creating a Leslie matrix from other forms of non-standard early
  age-classes can be accomplished by pre-processing \code{lx} and
  using the option \code{infant.class=FALSE}.

  The human sex ratio at birth (male births/female births) is
  remarkably close to SRB=1.05 across a wide range of populations and
  this is the default value for \code{SRB}.
 
  The resulting matrix has class "leslie.matrix". This class is
  not used extensively but will be in future development. 
 
 }

\value{
A k x k age-classified demographic projection matrix with class
"leslie.matrix". }

\references{ 

Keyfitz, N. 1977. Introduction to the mathematics of populations. 2nd
ed. Menlo Park: Addison-Wesley.

Preston, S.H., P. Heuveline, and F. Guillot. 2001. Demography:
Measuring and modeling population processes. Oxford: Blackwell. 

Caswell, H. 2001. Matrix population models: Construction, analysis,
and interpretation. 2nd ed. Sunderland, MA: Sinauer. 

 }


\seealso{  \code{\link{odiag}}, \code{\link{leslie.row1}}  }
\examples{

## Construct lx and mx values
data(goodman)
mlt <- with(goodman, life.table(x=age, nDx=mad.nDx, nKx=mad.nKx))
mx <- goodman$mad.bx/goodman$mad.nKx

## Now make the Leslie matrix
mad <- leslie.matrix(lx=mlt$nLx, mx=mx)

## Might as well do something with it too!
eigen.analysis(mad)
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
\keyword{ algebra }
