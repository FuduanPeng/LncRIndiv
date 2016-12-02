\name{LncRIndiv-package}
\alias{LncRIndiv-package}
\alias{LncRIndiv}
\docType{package}
\title{
Individual-level differential expression detecting method
}
\description{
To detect differentially expressed lncRNAs at the individual level.
}
\usage{
LncRIndiv(lncRNA, normal_expr, cancer_expr, freq)
}
\arguments{
  \item{lncRNA}{
  a (non-empty) vector of lncRNAs. The length of lncRNAs must be the same as the number of rows in the normal_expr and cancer_expr.
}
  \item{normal_expr}{
  a (non-empty) numeric matrix(or table) of data values.
}
  \item{cancer_expr}{
  a (non-empty) numeric matrix(or table) of data values.
}
\item{freq}{
  the criteria of identifying stable lncRNA pairs in normal(control) samples. The defult setting of freq is 0.95.
}
}
\details{
The fisher's exaxt test is used to test whether the stable lncRNA pair is significantly reversed in the cancer samples compared to normal samples. 
}
\value{
a data.frame containing five components is returned:
	\item{DE lncRNAs}{
		the lncRNAs that are differentially expressed in individual sample.
	}
	\item{paired lncRNAs}{
		the lncRNAs that are used to identify the DE lncRNAs.
	}
	\item{DE direction}{
		the dysregulation direction (up- or down-regulated) for each DE lncRNA.
	}
	\item{P-value}{
		the p-value of the Fisher's exaxt test.
	}
	\item{FDR}{
		the false discovery rate of P-value based on Benjamini & Hochbery.
	}
}
\author{
Maintainer: Fuduan Peng <pengfd@ems.hrbmu.edu.cn>, zjzhaozhangxiang@126.com and Yunyan Gu <guyunyan@ems.hrbmu.edu.cn>
}
\references{
Differential expression analysis at the individual level reveals a lncRNA prognostic signature for lung adenocarcinoma. Peng et. al.
}
\keyword{ package }
\seealso{
fisher.test: Fisher's exact test for testing the null of independence of rows and columns in a contingency table with fixed marginals.
}
\examples{
pair_info <- LncRIndiv(lncRNA,normal_expr,cancer_expr, freq)
}