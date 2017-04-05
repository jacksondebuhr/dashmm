#ifndef __TRACEUTILS_SQLITE_3_OUT_H__
#define __TRACEUTILS_SQLITE_3_OUT_H__


#include <string>
#include <vector>

#include <sqlite3.h>

#include "segment.h"
#include "trace.h"


namespace traceutils {


class SQLiteStatement {
 public:
  SQLiteStatement(sqlite3 *db, const std::string &sql);
  ~SQLiteStatement();

  // returns true if there are more steps that should be taken
  bool Step();

  // if clear is true, the bindings are cleared as well
  void Reset(bool clear = true);

  // Note that param indices start at 1
  int ParameterCount();
  void Bind(int index, double value);
  void Bind(int index, sqlite3_int64 value);
  void Bind(int index);  // bind null
  void Bind(int index, const std::string &value);

  // Get column info
  // Note that column indices start at 0
  //int ColumnCount();
  // How to overload...

 private:
  sqlite3 *db_;
  sqlite3_stmt *s_;
  std::string sql_;
};


class SQLiteWriter {
public:
  SQLiteWriter(const std::string &fname, const Trace &trace);
  ~SQLiteWriter();

  // this object will handle dealing with the statements and things.
  // This will present an interface that makes it easy to write out the
  // trace in the db format.
  void write();

private:
  // Likely a mess of private routines here
  void LocalityTable();
  void WorkerTable();
  void SegmenttypeTable();
  std::map<std::string, int> EventclassTable();
  std::map<std::string, int> EventtypeTable(
                                      const std::map<std::string, int> &cls);
  SegmentCollator EventTable(const std::map<std::string, int> &evttyp);
  void SegmentTable(const SegmentCollator &segs);

  std::string file_name_;
  sqlite3 *db_;
  const Trace &trace_;
};


} // traceutils


#endif // __TRACEUTILS_SQLITE_3_OUT_H__