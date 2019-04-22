package com.example.myapp;

import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import java.util.List;

import android.app.Activity;
import android.content.pm.ApplicationInfo;
import android.content.pm.PackageManager;
import android.os.Bundle;
import android.widget.ListView;

public class MainActivity extends Activity {

    private ListView appList;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        appList = (ListView) findViewById(R.id.appList);

        PackageManager pm = getPackageManager();
        List<ApplicationInfo> packages = pm
                .getInstalledApplications(PackageManager.GET_META_DATA);

        AppListAdapter appListAdapter = new AppListAdapter(this, pm, packages);
        appList.setAdapter(appListAdapter);
    }


}
/*public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
    }
}*/
