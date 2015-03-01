mpi=list(
rec=list(
pgi=list(
name='recruitment', nIVs=1, nPVs=6, nXIs=0, nPCs=1, 
pcs=list(
'1'=list(ids.PC=structure(c(1, 1, 1, 1, 1, 1, 1),names=c('YEAR_BLOCK', 'pLnR', 'pLnRCV', 'pLgtRX', 'pLnRa', 'pLnRb', 'pDevsLnR'),dim=c(7)), 
ids.Mod=structure(c(1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013),
    dimnames=list(index=c(1:64), type=c('YEAR_BLOCK')),dim=c(64, 1))),
NULL)
), 
pLnR=list(
  '1'=list(lower=0, upper=10, jitter='ON', initVal=4.3, finalVal=4.17397e-311, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnRCV=list(
  '1'=list(lower=-2, upper=2, jitter='OFF', initVal=-0.693147, finalVal=1.49167e-154, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLgtRX=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=0, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnRa=list(
  '1'=list(lower=1, upper=4, jitter='OFF', initVal=2.44235, finalVal=0, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnRb=list(
  '1'=list(lower=0, upper=4, jitter='OFF', initVal=1.38629, finalVal=0, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pDevsLnR=list(
  '1'=list(lower=-10, upper=10, jitter='OFF', initVal=0, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL), 
initVals=structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),names=c('1950', '1951', '1952', '1953', '1954', '1955', '1956', '1957', '1958', '1959', '1960', '1961', '1962', '1963', '1964', '1965', '1966', '1967', '1968', '1969', '1970', '1971', '1972', '1973', '1974', '1975', '1976', '1977', '1978', '1979', '1980', '1981', '1982', '1983', '1984', '1985', '1986', '1987', '1988', '1989', '1990', '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'),dim=c(64)), 
finalVals=structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),names=c('1950', '1951', '1952', '1953', '1954', '1955', '1956', '1957', '1958', '1959', '1960', '1961', '1962', '1963', '1964', '1965', '1966', '1967', '1968', '1969', '1970', '1971', '1972', '1973', '1974', '1975', '1976', '1977', '1978', '1979', '1980', '1981', '1982', '1983', '1984', '1985', '1986', '1987', '1988', '1989', '1990', '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'),dim=c(64))))

), 
nm=list(
pgi=list(
name='natural_mortality', nIVs=1, nPVs=5, nXIs=1, nPCs=1, 
pcs=list(
'1'=list(ids.PC=structure(c(1, 1, 0, 1, 1, 1, 0),names=c('YEAR_BLOCK', 'pLnM', 'pLnDMT', 'pLnDMX', 'pLnDMM', 'pLnDMXM', 'zScaling'),dim=c(7)), 
ids.Mod=structure(c(1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013),
    dimnames=list(index=c(1:64), type=c('YEAR_BLOCK')),dim=c(64, 1))),
NULL)
), 
pLnM=list(
  '1'=list(lower=-3, upper=2, jitter='OFF', initVal=-1.46968, finalVal=4.17397e-311, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDXT=list(
  '1'=list(lower=-2, upper=2, jitter='OFF', initVal=0, finalVal=1.49167e-154, phase=-4, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDMX=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=1.49167e-154, phase=-6, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDMM=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=0, phase=-6, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDMXM=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=0, phase=-6, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))

), 
grw=list(
pgi=list(
name='growth', nIVs=2, nPVs=3, nXIs=0, nPCs=2, 
pcs=list(
'1'=list(ids.PC=structure(c(1, 1, 1, 1, 1),names=c('YEAR_BLOCK', 'SEX', 'pLnGrA', 'pLnGrB', 'pLnGrBeta'),dim=c(5)), 
ids.Mod=structure(c(1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    dimnames=list(index=c(1:64), type=c('YEAR_BLOCK', 'SEX')),dim=c(64, 2))),
'2'=list(ids.PC=structure(c(1, 2, 2, 2, 1),names=c('YEAR_BLOCK', 'SEX', 'pLnGrA', 'pLnGrB', 'pLnGrBeta'),dim=c(5)), 
ids.Mod=structure(c(1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
    dimnames=list(index=c(1:64), type=c('YEAR_BLOCK', 'SEX')),dim=c(64, 2))),
NULL)
), 
pLnGrA=list(
  '1'=list(lower=-1.20397, upper=-0.510826, jitter='OFF', initVal=-0.798508, finalVal=4.1485e-311, phase=-7, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '2'=list(lower=-0.916291, upper=-0.356675, jitter='OFF', initVal=-0.597837, finalVal=1.49167e-154, phase=-7, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnGrB=list(
  '1'=list(lower=-0.356675, upper=0.182322, jitter='OFF', initVal=-0.0512933, finalVal=1.49167e-154, phase=-7, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '2'=list(lower=-0.510826, upper=0.182322, jitter='OFF', initVal=-0.105361, finalVal=0, phase=-7, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnGrBeta=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=-0.287682, finalVal=0, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))

), 
mat=list(
pgi=list(
name='maturity', nIVs=2, nPVs=1, nXIs=0, nPCs=2, 
pcs=list(
'1'=list(ids.PC=structure(c(1, 1, 1),names=c('YEAR_BLOCK', 'SEX', 'pvLgtPrMat'),dim=c(3)), 
ids.Mod=structure(c(1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    dimnames=list(index=c(1:64), type=c('YEAR_BLOCK', 'SEX')),dim=c(64, 2))),
'2'=list(ids.PC=structure(c(1, 2, 2),names=c('YEAR_BLOCK', 'SEX', 'pvLgtPrMat'),dim=c(3)), 
ids.Mod=structure(c(1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
    dimnames=list(index=c(1:64), type=c('YEAR_BLOCK', 'SEX')),dim=c(64, 2))),
NULL)
), 
pLgtPrMat=list(
  '1'=list(lower=-10, upper=10, jitter='OFF', initVal=0, finalVal=5.03947e-322, phase=-2, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL), 
initVals=structure(c(-7.25, -6.75, -6.25, -5.75, -5.25, -4.75, -4.25, -3.75, -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25),names=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32'),dim=c(32)), 
finalVals=structure(c(-7.25, -6.75, -6.25, -5.75, -5.25, -4.75, -4.25, -3.75, -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25),names=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32'),dim=c(32))),
  '2'=list(lower=-10, upper=10, jitter='OFF', initVal=0, finalVal=1.49167e-154, phase=-2, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL), 
initVals=structure(c(-5.25, -4.75, -4.25, -3.75, -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 9.99998),names=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32'),dim=c(32)), 
finalVals=structure(c(-5.25, -4.75, -4.25, -3.75, -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 9.99998),names=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32'),dim=c(32))))

), 
sel=list(
pgi=list(
name='selectivities', nIVs=1, nPVs=12, nXIs=2, nPCs=4, 
pcs=list(
'1'=list(ids.PC=structure(c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 182, 4),names=c('YEAR_BLOCK', 'pS1', 'pS2', 'pS3', 'pS4', 'pS5', 'pS6', 'pDevsS1', 'pDevsS2', 'pDevsS3', 'pDevsS4', 'pDevsS5', 'pDevsS6', 'fsZ', 'selFcn'),dim=c(15)), 
ids.Mod=structure(c(1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014),
    dimnames=list(index=c(1:40), type=c('YEAR_BLOCK')),dim=c(40, 1))),
'2'=list(ids.PC=structure(c(1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 182, 4),names=c('YEAR_BLOCK', 'pS1', 'pS2', 'pS3', 'pS4', 'pS5', 'pS6', 'pDevsS1', 'pDevsS2', 'pDevsS3', 'pDevsS4', 'pDevsS5', 'pDevsS6', 'fsZ', 'selFcn'),dim=c(15)), 
ids.Mod=structure(c(1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013),
    dimnames=list(index=c(1:49), type=c('YEAR_BLOCK')),dim=c(49, 1))),
'3'=list(ids.PC=structure(c(1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 182, 4),names=c('YEAR_BLOCK', 'pS1', 'pS2', 'pS3', 'pS4', 'pS5', 'pS6', 'pDevsS1', 'pDevsS2', 'pDevsS3', 'pDevsS4', 'pDevsS5', 'pDevsS6', 'fsZ', 'selFcn'),dim=c(15)), 
ids.Mod=structure(c(1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013),
    dimnames=list(index=c(1:49), type=c('YEAR_BLOCK')),dim=c(49, 1))),
'4'=list(ids.PC=structure(c(1, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 182, 4),names=c('YEAR_BLOCK', 'pS1', 'pS2', 'pS3', 'pS4', 'pS5', 'pS6', 'pDevsS1', 'pDevsS2', 'pDevsS3', 'pDevsS4', 'pDevsS5', 'pDevsS6', 'fsZ', 'selFcn'),dim=c(15)), 
ids.Mod=structure(c(1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013),
    dimnames=list(index=c(1:49), type=c('YEAR_BLOCK')),dim=c(49, 1))),
NULL)
), 
pS1=list(
  '1'=list(lower=1.6, upper=5.2, jitter='ON', initVal=4.09434, finalVal=3.26083e-322, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '2'=list(lower=4.4, upper=6, jitter='ON', initVal=4.94164, finalVal=1.49167e-154, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '3'=list(lower=3.5, upper=5.1, jitter='ON', initVal=4.38203, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '4'=list(lower=3.5, upper=5.1, jitter='ON', initVal=4.38203, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pS2=list(
  '1'=list(lower=1.6, upper=4, jitter='ON', initVal=3.2, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '2'=list(lower=1.6, upper=4, jitter='ON', initVal=2.3, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '3'=list(lower=1.6, upper=4, jitter='ON', initVal=3.2, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)),
  '4'=list(lower=1.6, upper=4, jitter='ON', initVal=3.2, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pS3=NULL, 
pS4=NULL, 
pS5=NULL, 
pS6=NULL, 
pDevsS1=NULL, 
pDevsS2=NULL, 
pDevsS3=NULL, 
pDevsS4=NULL, 
pDevsS5=NULL, 
pDevsS6=NULL
), 
fsh=list(
pgi=list(
name='fisheries', nIVs=3, nPVs=7, nXIs=3, nPCs=2, 
pcs=list(
'1'=list(ids.PC=structure(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 3, 2, 0),names=c('FISHERY', 'YEAR_BLOCK', 'SEX', 'pHM', 'pLnC', 'pLnDCT', 'pLnDCX', 'pLnDCM', 'pLnDCXM', 'pDevsLnC', 'idx.SelFcn', 'idx.RetFcn', 'useEffortRatio'),dim=c(13)), 
ids.Mod=structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    dimnames=list(index=c(1:49), type=c('FISHERY', 'YEAR_BLOCK', 'SEX')),dim=c(49, 3))),
'2'=list(ids.PC=structure(c(1, 1, 2, 1, 1, 0, 0, 0, 0, 1, 4, 0, 0),names=c('FISHERY', 'YEAR_BLOCK', 'SEX', 'pHM', 'pLnC', 'pLnDCT', 'pLnDCX', 'pLnDCM', 'pLnDCXM', 'pDevsLnC', 'idx.SelFcn', 'idx.RetFcn', 'useEffortRatio'),dim=c(13)), 
ids.Mod=structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
    dimnames=list(index=c(1:49), type=c('FISHERY', 'YEAR_BLOCK', 'SEX')),dim=c(49, 3))),
NULL)
), 
pLnC=list(
  '1'=list(lower=-5, upper=2, jitter='ON', initVal=-1.20397, finalVal=4.23975e-311, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDXT=list(
  '1'=list(lower=-5, upper=2, jitter='OFF', initVal=0, finalVal=1.49167e-154, phase=-2, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDCX=list(
  '1'=list(lower=-5, upper=2, jitter='OFF', initVal=0, finalVal=0, phase=-2, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDCM=list(
  '1'=list(lower=-5, upper=2, jitter='OFF', initVal=0, finalVal=0, phase=-2, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDCXM=list(
  '1'=list(lower=-5, upper=2, jitter='OFF', initVal=0, finalVal=0, phase=-2, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pDevsLnC=list(
  '1'=list(lower=-15, upper=15, jitter='OFF', initVal=0, finalVal=0, phase=1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL), 
initVals=structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),names=c('1965', '1966', '1967', '1968', '1969', '1970', '1971', '1972', '1973', '1974', '1975', '1976', '1977', '1978', '1979', '1980', '1981', '1982', '1983', '1984', '1985', '1986', '1987', '1988', '1989', '1990', '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'),dim=c(49)), 
finalVals=structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),names=c('1965', '1966', '1967', '1968', '1969', '1970', '1971', '1972', '1973', '1974', '1975', '1976', '1977', '1978', '1979', '1980', '1981', '1982', '1983', '1984', '1985', '1986', '1987', '1988', '1989', '1990', '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013'),dim=c(49))))

), 
srv=list(
pgi=list(
name='surveys', nIVs=3, nPVs=5, nXIs=1, nPCs=2, 
pcs=list(
'1'=list(ids.PC=structure(c(1, 1, 1, 1, 0, 0, 0, 0, 1),names=c('SURVEY', 'YEAR_BLOCK', 'SEX', 'pLnQ', 'pLnDQT', 'pLnDQX', 'pLnDQM', 'pLnDQXM', 'idx.SelFcn'),dim=c(9)), 
ids.Mod=structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    dimnames=list(index=c(1:40), type=c('SURVEY', 'YEAR_BLOCK', 'SEX')),dim=c(40, 3))),
'2'=list(ids.PC=structure(c(1, 1, 2, 1, 0, 1, 0, 0, 1),names=c('SURVEY', 'YEAR_BLOCK', 'SEX', 'pLnQ', 'pLnDQT', 'pLnDQX', 'pLnDQM', 'pLnDQXM', 'idx.SelFcn'),dim=c(9)), 
ids.Mod=structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
    dimnames=list(index=c(1:40), type=c('SURVEY', 'YEAR_BLOCK', 'SEX')),dim=c(40, 3))),
NULL)
), 
pLnQ=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=6.81811e-322, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDXT=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=1.49167e-154, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDQX=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=3.35965e-322, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDQM=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=1.49167e-154, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))
, 
pLnDQXM=list(
  '1'=list(lower=-1, upper=1, jitter='OFF', initVal=0, finalVal=0, phase=-1, resample='OFF', priorWgt=1, pdfType=list(type='none',params=NULL,consts=NULL)))

)
)